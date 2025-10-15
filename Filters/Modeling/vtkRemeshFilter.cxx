#include "vtkRemeshFilter.h"

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkFieldData.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkMath.h>
#include <vtkObjectFactory.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkTriangleFilter.h>

#include <vtkNew.h>

#include "vtk_eigen.h"
#include VTK_EIGEN(Core)
#include VTK_EIGEN(Dense)

#include <algorithm>
#include <cstdint>
#include <cmath>
#include <array>
#include <limits>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace
{

using Vector3d = Eigen::Vector3d;
using Vector3i = Eigen::Vector3i;

struct EdgeKey
{
  int A;
  int B;

  bool operator==(const EdgeKey& other) const noexcept
  {
    return this->A == other.A && this->B == other.B;
  }
};

struct EdgeKeyHash
{
  std::size_t operator()(const EdgeKey& key) const noexcept
  {
    return (static_cast<std::size_t>(static_cast<unsigned int>(key.A)) << 32) ^
      static_cast<unsigned int>(key.B);
  }
};

struct EdgeRecord
{
  int V0 = -1;
  int V1 = -1;
  std::vector<int> Faces;
};

struct MeshData
{
  std::vector<Vector3d> Vertices;
  std::vector<Vector3i> Faces;
};

struct VertexAttributes
{
  std::vector<char> IsFeature;
  std::vector<double> High;
  std::vector<double> Low;
};

constexpr double DegenerateAreaThreshold = 1e-18;

EdgeKey MakeEdgeKey(int a, int b)
{
  if (a > b)
  {
    std::swap(a, b);
  }
  return EdgeKey{ a, b };
}

struct FaceKey
{
  int V[3];

  bool operator==(const FaceKey& other) const noexcept
  {
    return this->V[0] == other.V[0] && this->V[1] == other.V[1] && this->V[2] == other.V[2];
  }
};

struct FaceKeyHash
{
  std::size_t operator()(const FaceKey& key) const noexcept
  {
    std::size_t h0 = static_cast<std::size_t>(static_cast<unsigned int>(key.V[0]));
    std::size_t h1 = static_cast<std::size_t>(static_cast<unsigned int>(key.V[1]));
    std::size_t h2 = static_cast<std::size_t>(static_cast<unsigned int>(key.V[2]));
    return ((h0 * 73856093u) ^ (h1 * 19349663u) ^ (h2 * 83492791u));
  }
};

bool RemoveDegenerateFaces(MeshData& mesh)
{
  bool removed = false;
  std::vector<Vector3i> newFaces;
  newFaces.reserve(mesh.Faces.size());
  std::unordered_set<FaceKey, FaceKeyHash> uniqueFaces;
  uniqueFaces.reserve(mesh.Faces.size());

  for (const auto& face : mesh.Faces)
  {
    int a = face[0];
    int b = face[1];
    int c = face[2];
    if (a < 0 || b < 0 || c < 0 || a == b || b == c || c == a)
    {
      removed = true;
      continue;
    }

    const Vector3d& p0 = mesh.Vertices[a];
    const Vector3d& p1 = mesh.Vertices[b];
    const Vector3d& p2 = mesh.Vertices[c];
    Vector3d normal = (p1 - p0).cross(p2 - p0);
    if (normal.squaredNorm() <= DegenerateAreaThreshold)
    {
      removed = true;
      continue;
    }

    std::array<int, 3> sorted = { a, b, c };
    std::sort(sorted.begin(), sorted.end());
    FaceKey key{ { sorted[0], sorted[1], sorted[2] } };
    if (!uniqueFaces.insert(key).second)
    {
      removed = true;
      continue;
    }

    newFaces.push_back(face);
  }

  if (removed)
  {
    mesh.Faces.swap(newFaces);
  }

  return removed;
}

std::unordered_map<EdgeKey, EdgeRecord, EdgeKeyHash> BuildEdgeRecords(const MeshData& mesh)
{
  std::unordered_map<EdgeKey, EdgeRecord, EdgeKeyHash> edges;
  edges.reserve(mesh.Faces.size() * 3);

  for (std::size_t faceIndex = 0; faceIndex < mesh.Faces.size(); ++faceIndex)
  {
    const Vector3i& face = mesh.Faces[faceIndex];
    auto addEdge = [&](int u, int v) {
      EdgeKey key = MakeEdgeKey(u, v);
      EdgeRecord& record = edges[key];
      record.V0 = key.A;
      record.V1 = key.B;
      record.Faces.push_back(static_cast<int>(faceIndex));
    };

    addEdge(face[0], face[1]);
    addEdge(face[1], face[2]);
    addEdge(face[2], face[0]);
  }

  return edges;
}

std::vector<std::unordered_set<int>> BuildVertexNeighborSets(const MeshData& mesh)
{
  std::vector<std::unordered_set<int>> neighbors(mesh.Vertices.size());
  for (const auto& face : mesh.Faces)
  {
    int a = face[0];
    int b = face[1];
    int c = face[2];
    neighbors[a].insert(b);
    neighbors[a].insert(c);
    neighbors[b].insert(a);
    neighbors[b].insert(c);
    neighbors[c].insert(a);
    neighbors[c].insert(b);
  }
  return neighbors;
}

void CompactMesh(MeshData& mesh, VertexAttributes& attr, const std::vector<char>& removed)
{
  if (mesh.Vertices.empty())
  {
    return;
  }

  std::vector<int> mapping(mesh.Vertices.size(), -1);
  MeshData compacted;
  compacted.Vertices.reserve(mesh.Vertices.size());
  compacted.Faces.reserve(mesh.Faces.size());

  VertexAttributes newAttr;
  newAttr.IsFeature.reserve(attr.IsFeature.size());
  newAttr.High.reserve(attr.High.size());
  newAttr.Low.reserve(attr.Low.size());

  for (std::size_t i = 0; i < mesh.Vertices.size(); ++i)
  {
    if (!removed[i])
    {
      mapping[i] = static_cast<int>(compacted.Vertices.size());
      compacted.Vertices.push_back(mesh.Vertices[i]);
      newAttr.IsFeature.push_back(attr.IsFeature[i]);
      newAttr.High.push_back(attr.High[i]);
      newAttr.Low.push_back(attr.Low[i]);
    }
  }

  for (const auto& face : mesh.Faces)
  {
    int a = mapping[face[0]];
    int b = mapping[face[1]];
    int c = mapping[face[2]];
    if (a < 0 || b < 0 || c < 0 || a == b || b == c || c == a)
    {
      continue;
    }
    compacted.Faces.push_back(Vector3i(a, b, c));
  }

  mesh = std::move(compacted);
  attr = std::move(newAttr);
}

std::vector<char> DetectFeatureVertices(const MeshData& mesh, double angleRadians)
{
  std::vector<char> features(mesh.Vertices.size(), 0);
  if (mesh.Faces.empty())
  {
    return features;
  }

  auto edges = BuildEdgeRecords(mesh);
  double cosThreshold = std::cos(angleRadians);

  auto faceNormal = [&](const Vector3i& tri) {
    const Vector3d& p0 = mesh.Vertices[tri[0]];
    const Vector3d& p1 = mesh.Vertices[tri[1]];
    const Vector3d& p2 = mesh.Vertices[tri[2]];
    Vector3d n = (p1 - p0).cross(p2 - p0);
    double len = n.norm();
    if (len > 0.0)
    {
      n /= len;
    }
    return n;
  };

  for (const auto& entry : edges)
  {
    const EdgeRecord& record = entry.second;
    if (record.Faces.size() != 2)
    {
      continue;
    }

    const Vector3i& f0 = mesh.Faces[record.Faces[0]];
    const Vector3i& f1 = mesh.Faces[record.Faces[1]];
    Vector3d n0 = faceNormal(f0);
    Vector3d n1 = faceNormal(f1);
    double dot = std::clamp(n0.dot(n1), -1.0, 1.0);
    if (dot <= cosThreshold)
    {
      features[record.V0] = 1;
      features[record.V1] = 1;
    }
  }

  return features;
}

bool TrySplitEdge(MeshData& mesh, VertexAttributes& attr)
{
  auto edges = BuildEdgeRecords(mesh);
  for (const auto& entry : edges)
  {
    int v0 = entry.second.V0;
    int v1 = entry.second.V1;
    if (v0 < 0 || v1 < 0)
    {
      continue;
    }
    if (attr.IsFeature[v0] || attr.IsFeature[v1])
    {
      continue;
    }

    double length = (mesh.Vertices[v0] - mesh.Vertices[v1]).norm();
    double threshold = 0.5 * (attr.High[v0] + attr.High[v1]);
    if (length <= threshold)
    {
      continue;
    }

    int newIndex = static_cast<int>(mesh.Vertices.size());
    mesh.Vertices.push_back(0.5 * (mesh.Vertices[v0] + mesh.Vertices[v1]));
    attr.IsFeature.push_back(0);
    attr.High.push_back(0.5 * (attr.High[v0] + attr.High[v1]));
    attr.Low.push_back(0.5 * (attr.Low[v0] + attr.Low[v1]));

    std::vector<Vector3i> newFaces;
    for (int faceIndex : entry.second.Faces)
    {
      if (faceIndex < 0 || faceIndex >= static_cast<int>(mesh.Faces.size()))
      {
        continue;
      }

      Vector3i& face = mesh.Faces[faceIndex];
      Vector3i originalFace = face;
      int pos0 = -1;
      int pos1 = -1;
      for (int k = 0; k < 3; ++k)
      {
        if (face[k] == v0)
        {
          pos0 = k;
        }
        else if (face[k] == v1)
        {
          pos1 = k;
        }
      }

      if (pos0 < 0 || pos1 < 0)
      {
        continue;
      }

      int pos2 = 3 - pos0 - pos1;
      int v2 = face[pos2];

      Vector3i faceA;
      Vector3i faceB;

      if ((pos0 + 1) % 3 == pos1)
      {
        faceA = Vector3i(v0, newIndex, v2);
        faceB = Vector3i(newIndex, v1, v2);
      }
      else
      {
        faceA = Vector3i(v1, newIndex, v2);
        faceB = Vector3i(newIndex, v0, v2);
      }

      Vector3d originalNormal =
        (mesh.Vertices[originalFace[1]] - mesh.Vertices[originalFace[0]]).cross(
          mesh.Vertices[originalFace[2]] - mesh.Vertices[originalFace[0]]);

      Vector3d normalA = (mesh.Vertices[faceA[1]] - mesh.Vertices[faceA[0]]).cross(
        mesh.Vertices[faceA[2]] - mesh.Vertices[faceA[0]]);
      if (normalA.dot(originalNormal) < 0.0)
      {
        faceA = Vector3i(faceA[0], faceA[2], faceA[1]);
      }

      Vector3d normalB = (mesh.Vertices[faceB[1]] - mesh.Vertices[faceB[0]]).cross(
        mesh.Vertices[faceB[2]] - mesh.Vertices[faceB[0]]);
      if (normalB.dot(originalNormal) < 0.0)
      {
        faceB = Vector3i(faceB[0], faceB[2], faceB[1]);
      }

      face = faceA;
      newFaces.push_back(faceB);
    }

    for (const auto& f : newFaces)
    {
      mesh.Faces.push_back(f);
    }

    RemoveDegenerateFaces(mesh);
    return true;
  }

  return false;
}

bool SplitEdgesUntilBound(MeshData& mesh, VertexAttributes& attr)
{
  bool changed = false;
  while (TrySplitEdge(mesh, attr))
  {
    changed = true;
  }
  return changed;
}

bool TryCollapseEdge(MeshData& mesh, VertexAttributes& attr)
{
  auto edges = BuildEdgeRecords(mesh);
  std::vector<char> boundary(mesh.Vertices.size(), 0);
  for (const auto& entry : edges)
  {
    if (entry.second.Faces.size() == 1)
    {
      boundary[entry.second.V0] = 1;
      boundary[entry.second.V1] = 1;
    }
  }

  for (const auto& entry : edges)
  {
    int v0 = entry.second.V0;
    int v1 = entry.second.V1;
    if (v0 < 0 || v1 < 0)
    {
      continue;
    }
    if (attr.IsFeature[v0] || attr.IsFeature[v1])
    {
      continue;
    }
    if (boundary[v0] != boundary[v1])
    {
      continue;
    }

    double length = (mesh.Vertices[v0] - mesh.Vertices[v1]).norm();
    double threshold = 0.5 * (attr.Low[v0] + attr.Low[v1]);
    if (length >= threshold)
    {
      continue;
    }

    int keep = std::min(v0, v1);
    int remove = std::max(v0, v1);

    mesh.Vertices[keep] = 0.5 * (mesh.Vertices[keep] + mesh.Vertices[remove]);
    attr.High[keep] = 0.5 * (attr.High[keep] + attr.High[remove]);
    attr.Low[keep] = 0.5 * (attr.Low[keep] + attr.Low[remove]);
    attr.IsFeature[keep] = attr.IsFeature[keep] || attr.IsFeature[remove];

    std::vector<char> removed(mesh.Vertices.size(), 0);
    removed[remove] = 1;

    for (auto& face : mesh.Faces)
    {
      for (int k = 0; k < 3; ++k)
      {
        if (face[k] == remove)
        {
          face[k] = keep;
        }
      }
    }

    RemoveDegenerateFaces(mesh);
    CompactMesh(mesh, attr, removed);
    RemoveDegenerateFaces(mesh);
    return true;
  }

  return false;
}

bool CollapseEdges(MeshData& mesh, VertexAttributes& attr)
{
  bool changed = false;
  while (TryCollapseEdge(mesh, attr))
  {
    changed = true;
  }
  return changed;
}

bool TryFlipEdge(MeshData& mesh, const VertexAttributes& attr)
{
  auto edges = BuildEdgeRecords(mesh);
  auto neighbors = BuildVertexNeighborSets(mesh);

  std::vector<char> boundary(mesh.Vertices.size(), 0);
  for (const auto& entry : edges)
  {
    if (entry.second.Faces.size() == 1)
    {
      boundary[entry.second.V0] = 1;
      boundary[entry.second.V1] = 1;
    }
  }

  auto deviation = [&](int v, int valence) {
    int target = boundary[v] ? 4 : 6;
    double diff = static_cast<double>(valence - target);
    return diff * diff;
  };

  for (const auto& entry : edges)
  {
    const EdgeRecord& record = entry.second;
    if (record.Faces.size() != 2)
    {
      continue;
    }

    int v0 = record.V0;
    int v1 = record.V1;
    if (attr.IsFeature[v0] || attr.IsFeature[v1])
    {
      continue;
    }

    const Vector3i& face0 = mesh.Faces[record.Faces[0]];
    const Vector3i& face1 = mesh.Faces[record.Faces[1]];
    int c = -1;
    int d = -1;
    for (int k = 0; k < 3; ++k)
    {
      if (face0[k] != v0 && face0[k] != v1)
      {
        c = face0[k];
      }
      if (face1[k] != v0 && face1[k] != v1)
      {
        d = face1[k];
      }
    }

    if (c < 0 || d < 0 || c == d)
    {
      continue;
    }
    if (attr.IsFeature[c] || attr.IsFeature[d])
    {
      continue;
    }

    EdgeKey diagonal = MakeEdgeKey(c, d);
    auto diagonalIt = edges.find(diagonal);
    if (diagonalIt != edges.end() && diagonalIt->second.Faces.size() >= 2)
    {
      continue;
    }

    int valenceV0 = static_cast<int>(neighbors[v0].size());
    int valenceV1 = static_cast<int>(neighbors[v1].size());
    int valenceC = static_cast<int>(neighbors[c].size());
    int valenceD = static_cast<int>(neighbors[d].size());

    double scoreBefore = deviation(v0, valenceV0) + deviation(v1, valenceV1) +
      deviation(c, valenceC) + deviation(d, valenceD);

    double scoreAfter = deviation(v0, valenceV0 - 1) + deviation(v1, valenceV1 - 1) +
      deviation(c, valenceC + (neighbors[c].count(d) ? 0 : 1)) +
      deviation(d, valenceD + (neighbors[d].count(c) ? 0 : 1));

    if (scoreAfter >= scoreBefore - 1e-8)
    {
      continue;
    }

    Vector3d n0 = (mesh.Vertices[face0[1]] - mesh.Vertices[face0[0]]).cross(
      mesh.Vertices[face0[2]] - mesh.Vertices[face0[0]]);
    Vector3d n1 = (mesh.Vertices[face1[1]] - mesh.Vertices[face1[0]]).cross(
      mesh.Vertices[face1[2]] - mesh.Vertices[face1[0]]);

    Vector3i newFace0(v0, d, c);
    Vector3i newFace1(v1, c, d);

    Vector3d nf0 = (mesh.Vertices[newFace0[1]] - mesh.Vertices[newFace0[0]]).cross(
      mesh.Vertices[newFace0[2]] - mesh.Vertices[newFace0[0]]);
    if (nf0.dot(n0) < 0.0)
    {
      newFace0 = Vector3i(v0, c, d);
    }

    Vector3d nf1 = (mesh.Vertices[newFace1[1]] - mesh.Vertices[newFace1[0]]).cross(
      mesh.Vertices[newFace1[2]] - mesh.Vertices[newFace1[0]]);
    if (nf1.dot(n1) < 0.0)
    {
      newFace1 = Vector3i(v1, d, c);
    }

    mesh.Faces[record.Faces[0]] = newFace0;
    mesh.Faces[record.Faces[1]] = newFace1;
    RemoveDegenerateFaces(mesh);
    return true;
  }

  return false;
}

bool EqualizeValences(MeshData& mesh, const VertexAttributes& attr)
{
  bool changed = false;
  while (TryFlipEdge(mesh, attr))
  {
    changed = true;
  }
  return changed;
}

void TangentialRelaxation(MeshData& mesh, const VertexAttributes& attr)
{
  const std::size_t n = mesh.Vertices.size();
  if (n == 0)
  {
    return;
  }

  auto neighbors = BuildVertexNeighborSets(mesh);
  std::vector<Vector3d> normals(n, Vector3d::Zero());

  for (const auto& face : mesh.Faces)
  {
    const Vector3d& p0 = mesh.Vertices[face[0]];
    const Vector3d& p1 = mesh.Vertices[face[1]];
    const Vector3d& p2 = mesh.Vertices[face[2]];
    Vector3d nrm = (p1 - p0).cross(p2 - p0);
    for (int i = 0; i < 3; ++i)
    {
      normals[face[i]] += nrm;
    }
  }

  for (auto& normal : normals)
  {
    double len = normal.norm();
    if (len > 0.0)
    {
      normal /= len;
    }
    else
    {
      normal = Vector3d(0.0, 0.0, 1.0);
    }
  }

  std::vector<Vector3d> updated = mesh.Vertices;
  for (std::size_t i = 0; i < n; ++i)
  {
    if (attr.IsFeature[i] || neighbors[i].empty())
    {
      continue;
    }

    Vector3d centroid = Vector3d::Zero();
    for (int neighbor : neighbors[i])
    {
      centroid += mesh.Vertices[static_cast<std::size_t>(neighbor)];
    }
    centroid /= static_cast<double>(neighbors[i].size());

    Vector3d displacement = centroid - mesh.Vertices[i];
    const Vector3d& normal = normals[i];
    displacement -= displacement.dot(normal) * normal;
    updated[i] = mesh.Vertices[i] + displacement;
  }

  mesh.Vertices.swap(updated);
}

Vector3d ClosestPointOnTriangle(const Vector3d& p, const Vector3d& a, const Vector3d& b, const Vector3d& c)
{
  const Vector3d ab = b - a;
  const Vector3d ac = c - a;
  const Vector3d ap = p - a;
  double d1 = ab.dot(ap);
  double d2 = ac.dot(ap);
  if (d1 <= 0.0 && d2 <= 0.0)
  {
    return a;
  }

  const Vector3d bp = p - b;
  double d3 = ab.dot(bp);
  double d4 = ac.dot(bp);
  if (d3 >= 0.0 && d4 <= d3)
  {
    return b;
  }

  double vc = d1 * d4 - d3 * d2;
  if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0)
  {
    double v = d1 / (d1 - d3);
    return a + v * ab;
  }

  const Vector3d cp = p - c;
  double d5 = ab.dot(cp);
  double d6 = ac.dot(cp);
  if (d6 >= 0.0 && d5 <= d6)
  {
    return c;
  }

  double vb = d5 * d2 - d1 * d6;
  if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0)
  {
    double w = d2 / (d2 - d6);
    return a + w * ac;
  }

  double va = d3 * d6 - d5 * d4;
  if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0)
  {
    double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
    return b + w * (c - b);
  }

  double denom = 1.0 / (va + vb + vc);
  double v = vb * denom;
  double w = vc * denom;
  return a + ab * v + ac * w;
}

void ProjectToSurface(MeshData& mesh, const VertexAttributes& attr, const MeshData& original)
{
  if (original.Faces.empty())
  {
    return;
  }

  for (std::size_t i = 0; i < mesh.Vertices.size(); ++i)
  {
    if (attr.IsFeature[i])
    {
      continue;
    }

    const Vector3d& p = mesh.Vertices[i];
    double bestDistance2 = std::numeric_limits<double>::max();
    Vector3d bestPoint = p;

    for (const auto& face : original.Faces)
    {
      const Vector3d& a = original.Vertices[face[0]];
      const Vector3d& b = original.Vertices[face[1]];
      const Vector3d& c = original.Vertices[face[2]];
      Vector3d candidate = ClosestPointOnTriangle(p, a, b, c);
      double d2 = (candidate - p).squaredNorm();
      if (d2 < bestDistance2)
      {
        bestDistance2 = d2;
        bestPoint = candidate;
      }
    }

    mesh.Vertices[i] = bestPoint;
  }
}

void RemeshBotsch(MeshData& mesh, const MeshData& original, double targetEdgeLength,
  double featureAngleDegrees, int iterations, bool projectToSurface)
{
  if (mesh.Vertices.empty() || mesh.Faces.empty())
  {
    return;
  }

  VertexAttributes attr;
  attr.IsFeature.assign(mesh.Vertices.size(), 0);
  attr.High.assign(mesh.Vertices.size(), 1.4 * targetEdgeLength);
  attr.Low.assign(mesh.Vertices.size(), 0.7 * targetEdgeLength);

  if (featureAngleDegrees > 0.0)
  {
    double angleRadians = vtkMath::RadiansFromDegrees(featureAngleDegrees);
    std::vector<char> features = DetectFeatureVertices(mesh, angleRadians);
    for (std::size_t i = 0; i < features.size() && i < attr.IsFeature.size(); ++i)
    {
      attr.IsFeature[i] = features[i];
    }
  }

  iterations = std::max(iterations, 1);
  for (int iter = 0; iter < iterations; ++iter)
  {
    SplitEdgesUntilBound(mesh, attr);
    CollapseEdges(mesh, attr);
    EqualizeValences(mesh, attr);
    TangentialRelaxation(mesh, attr);
    RemoveDegenerateFaces(mesh);
    if (projectToSurface)
    {
      ProjectToSurface(mesh, attr, original);
    }
  }

  RemoveDegenerateFaces(mesh);
}

MeshData MeshFromPolyData(vtkPolyData* poly)
{
  MeshData mesh;
  if (!poly)
  {
    return mesh;
  }

  vtkPoints* points = poly->GetPoints();
  vtkCellArray* polys = poly->GetPolys();
  if (!points || !polys)
  {
    return mesh;
  }

  vtkIdType numberOfPoints = points->GetNumberOfPoints();
  mesh.Vertices.resize(static_cast<std::size_t>(numberOfPoints));
  double p[3];
  for (vtkIdType i = 0; i < numberOfPoints; ++i)
  {
    points->GetPoint(i, p);
    mesh.Vertices[static_cast<std::size_t>(i)] = Vector3d(p[0], p[1], p[2]);
  }

  polys->InitTraversal();
  vtkIdType npts = 0;
  const vtkIdType* pts = nullptr;
  while (polys->GetNextCell(npts, pts))
  {
    if (npts != 3)
    {
      continue;
    }
    mesh.Faces.push_back(Vector3i(static_cast<int>(pts[0]), static_cast<int>(pts[1]),
      static_cast<int>(pts[2])));
  }

  RemoveDegenerateFaces(mesh);
  return mesh;
}

void CopyMeshToPolyData(const MeshData& mesh, vtkPolyData* output)
{
  if (!output)
  {
    return;
  }

  vtkNew<vtkPoints> points;
  points->SetNumberOfPoints(static_cast<vtkIdType>(mesh.Vertices.size()));
  for (vtkIdType i = 0; i < static_cast<vtkIdType>(mesh.Vertices.size()); ++i)
  {
    const Vector3d& p = mesh.Vertices[static_cast<std::size_t>(i)];
    points->SetPoint(i, p.data());
  }

  vtkNew<vtkCellArray> polys;
  for (const auto& face : mesh.Faces)
  {
    vtkIdType tri[3] = { face[0], face[1], face[2] };
    polys->InsertNextCell(3, tri);
  }

  output->SetPoints(points);
  output->SetPolys(polys);
}

} // end anonymous namespace

VTK_ABI_NAMESPACE_BEGIN
vtkStandardNewMacro(vtkRemeshFilter);

vtkRemeshFilter::vtkRemeshFilter()
  : TargetEdgeLength(0.0)
  , FeatureAngle(30.0)
  , Iterations(5)
  , ProjectToSurface(true)
{
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

void vtkRemeshFilter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
  os << indent << "TargetEdgeLength: " << this->TargetEdgeLength << "\n";
  os << indent << "FeatureAngle: " << this->FeatureAngle << "\n";
  os << indent << "Iterations: " << this->Iterations << "\n";
  os << indent << "ProjectToSurface: " << (this->ProjectToSurface ? "On" : "Off") << "\n";
}

int vtkRemeshFilter::FillInputPortInformation(int port, vtkInformation* info)
{
  if (port == 0)
  {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
    return 1;
  }
  return 0;
}

int vtkRemeshFilter::RequestData(
  vtkInformation*, vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
  vtkPolyData* input = vtkPolyData::GetData(inputVector[0], 0);
  vtkPolyData* output = vtkPolyData::GetData(outputVector, 0);

  if (!input || !output)
  {
    return 0;
  }

  vtkNew<vtkTriangleFilter> triangulator;
  triangulator->SetInputData(input);
  triangulator->PassVertsOff();
  triangulator->PassLinesOff();
  triangulator->Update();

  vtkPolyData* workMesh = triangulator->GetOutput();
  MeshData mesh = MeshFromPolyData(workMesh);
  if (mesh.Vertices.empty() || mesh.Faces.empty())
  {
    vtkWarningMacro(<< "Remeshing requires a triangulated surface. Returning input.");
    output->ShallowCopy(input);
    return 1;
  }

  double edgeLength = this->TargetEdgeLength;
  if (edgeLength <= 0.0)
  {
    edgeLength = this->ComputeDefaultEdgeLength(workMesh);
  }

  if (edgeLength <= 0.0)
  {
    vtkWarningMacro(<< "Unable to determine target edge length. Returning input.");
    output->ShallowCopy(input);
    return 1;
  }

  MeshData original = mesh;
  RemeshBotsch(mesh, original, edgeLength, this->FeatureAngle, this->Iterations,
    this->ProjectToSurface);

  if (mesh.Vertices.empty() || mesh.Faces.empty())
  {
    vtkWarningMacro(<< "Remeshing produced an empty mesh. Returning input.");
    output->ShallowCopy(input);
    return 1;
  }

  output->Initialize();
  CopyMeshToPolyData(mesh, output);
  output->GetPointData()->Initialize();
  output->GetCellData()->Initialize();
  if (input->GetFieldData())
  {
    output->GetFieldData()->ShallowCopy(input->GetFieldData());
  }

  return 1;
}

double vtkRemeshFilter::ComputeDefaultEdgeLength(vtkPolyData* mesh) const
{
  if (!mesh || !mesh->GetPoints())
  {
    return 0.0;
  }

  vtkPoints* points = mesh->GetPoints();
  vtkCellArray* polys = mesh->GetPolys();
  if (!polys)
  {
    return 0.0;
  }

  auto encodeEdge = [](vtkIdType a, vtkIdType b) -> std::uint64_t {
    if (a > b)
    {
      std::swap(a, b);
    }
    return (static_cast<std::uint64_t>(a) << 32) | static_cast<std::uint64_t>(b);
  };

  std::unordered_set<std::uint64_t> uniqueEdges;
  uniqueEdges.reserve(static_cast<std::size_t>(polys->GetNumberOfCells()) * 3);

  polys->InitTraversal();
  vtkIdType npts = 0;
  const vtkIdType* pts = nullptr;
  double total = 0.0;
  vtkIdType edgeCount = 0;

  double p0[3];
  double p1[3];
  while (polys->GetNextCell(npts, pts))
  {
    if (npts < 2)
    {
      continue;
    }

    for (vtkIdType i = 0; i < npts; ++i)
    {
      vtkIdType id0 = pts[i];
      vtkIdType id1 = pts[(i + 1) % npts];
      const std::uint64_t key = encodeEdge(id0, id1);
      if (uniqueEdges.insert(key).second)
      {
        points->GetPoint(id0, p0);
        points->GetPoint(id1, p1);
        double dx = p0[0] - p1[0];
        double dy = p0[1] - p1[1];
        double dz = p0[2] - p1[2];
        total += std::sqrt(dx * dx + dy * dy + dz * dz);
        edgeCount++;
      }
    }
  }

  return edgeCount > 0 ? total / static_cast<double>(edgeCount) : 0.0;
}

VTK_ABI_NAMESPACE_END
