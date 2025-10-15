#include "MeshIO.h"
#include "Face.h"
#include "HalfEdge.h"
#include "Mesh.h"
#include "Vertex.h"

#include <vtkCellArray.h>
#include <vtkNew.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace vtkBotschKobbeltRemeshing
{
namespace meshio_detail
{
class Index
{
public:
  Index() = default;

  Index(int v, int vt, int vn)
    : position(v)
    , uv(vt)
    , normal(vn)
  {
  }

  bool operator<(const Index& other) const
  {
    if (position < other.position)
    {
      return true;
    }
    if (position > other.position)
    {
      return false;
    }
    if (uv < other.uv)
    {
      return true;
    }
    if (uv > other.uv)
    {
      return false;
    }
    if (normal < other.normal)
    {
      return true;
    }
    if (normal > other.normal)
    {
      return false;
    }

    return false;
  }

  int position{ -1 };
  int uv{ -1 };
  int normal{ -1 };
};

struct EdgeKey
{
  int A{ -1 };
  int B{ -1 };

  bool operator==(const EdgeKey& other) const noexcept
  {
    return this->A == other.A && this->B == other.B;
  }
};

EdgeKey MakeEdgeKey(int a, int b)
{
  if (a > b)
  {
    std::swap(a, b);
  }
  EdgeKey key;
  key.A = a;
  key.B = b;
  return key;
}

struct EdgeKeyHash
{
  std::size_t operator()(const EdgeKey& key) const noexcept
  {
    return (static_cast<std::size_t>(static_cast<unsigned int>(key.A)) << 32) ^
      static_cast<unsigned int>(key.B);
  }
};

Index parseFaceIndex(const std::string& token)
{
  std::stringstream in(token);
  std::string indexString;
  int indices[3] = { -1, -1, -1 };

  int i = 0;
  while (getline(in, indexString, '/'))
  {
    if (indexString != "\\")
    {
      std::stringstream ss(indexString);
      ss >> indices[i++];
    }
  }

  return Index(indices[0] - 1, indices[1] - 1, indices[2] - 1);
}

std::string stringRep(const Eigen::Vector3d& v)
{
  return std::to_string(v.x()) + " " + std::to_string(v.y()) + " " + std::to_string(v.z());
}

} // namespace meshio_detail

using meshio_detail::Index;
using meshio_detail::EdgeKey;
using meshio_detail::EdgeKeyHash;
using meshio_detail::MakeEdgeKey;
using meshio_detail::parseFaceIndex;
using meshio_detail::stringRep;

class MeshData
{
public:
  std::vector<Eigen::Vector3d> positions;
  std::vector<Eigen::Vector3d> uvs;
  std::vector<Eigen::Vector3d> normals;
  std::vector<std::vector<Index>> indices;
};

extern std::vector<HalfEdge> isolated;

void MeshIO::preallocateMeshElements(const MeshData& data, Mesh& mesh)
{
  std::unordered_set<EdgeKey, EdgeKeyHash> edges;
  edges.reserve(data.indices.size() * 3);
  for (const auto& f : data.indices)
  {
    for (std::size_t I = 0; I < f.size(); I++)
    {
      std::size_t J = (I + 1) % f.size();
      int i = f[I].position;
      int j = f[J].position;

      edges.insert(MakeEdgeKey(i, j));
    }
  }

  std::size_t nV = data.positions.size();
  std::size_t nE = edges.size();
  std::size_t nF = data.indices.size();
  std::size_t nHE = 2 * nE;
  std::size_t chi = nV - nE + nF;
  int nB = std::max(0, 2 - static_cast<int>(chi));

  mesh.halfEdges.clear();
  mesh.vertices.clear();
  mesh.edges.clear();
  mesh.faces.clear();
  mesh.boundaries.clear();

  mesh.halfEdges.reserve(3 * nHE);
  mesh.vertices.reserve(3 * nV);
  mesh.edges.reserve(3 * nE);
  mesh.faces.reserve(3 * (nF + nB));
}

void MeshIO::indexElements(Mesh& mesh)
{
  int index = 0;
  for (VertexIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++)
  {
    v->index = index;
    index++;
  }

  index = 0;
  for (EdgeIter e = mesh.edges.begin(); e != mesh.edges.end(); e++)
  {
    e->index = index;
    index++;
  }

  index = 0;
  for (HalfEdgeIter he = mesh.halfEdges.begin(); he != mesh.halfEdges.end(); he++)
  {
    he->index = index;
    index++;
  }

  index = 0;
  for (FaceIter f = mesh.faces.begin(); f != mesh.faces.end(); f++)
  {
    f->index = index;
    index++;
  }
}

void MeshIO::checkIsolatedVertices(const Mesh& mesh)
{
  for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++)
  {
    if (v->isIsolated())
    {
      std::cerr << "Warning: vertex " << v->index
                << " is isolated (not contained in any face)." << std::endl;
    }
  }
}

void MeshIO::checkNonManifoldVertices(const Mesh& mesh)
{
  std::unordered_map<std::string, int> vertexFaceMap;

  for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++)
  {
    HalfEdgeCIter he = f->he;
    do
    {
  vertexFaceMap[stringRep(he->vertex->position)]++;
      he = he->next;

    } while (he != f->he);
  }

  for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++)
  {
    int valence = 0;
    HalfEdgeCIter he = v->he;
    do
    {
      valence++;
      he = he->flip->next;

    } while (he != v->he);

    if (vertexFaceMap[stringRep(v->position)] != valence)
    {
      std::cerr << "Warning: vertex " << v->index << " is nonmanifold." << std::endl;
    }
  }
}

bool MeshIO::buildMesh(const MeshData& data, Mesh& mesh)
{
  const std::size_t estimatedHalfEdgeCount = data.indices.size() * 3;
  std::unordered_map<EdgeKey, int, EdgeKeyHash> edgeCount;
  std::unordered_map<EdgeKey, int, EdgeKeyHash> existingHalfEdges;
  std::vector<VertexIter> indexToVertex(data.positions.size());
  std::vector<char> halfEdgeHasFlip;
  edgeCount.reserve(estimatedHalfEdgeCount);
  existingHalfEdges.reserve(estimatedHalfEdgeCount);
  halfEdgeHasFlip.reserve(estimatedHalfEdgeCount * 2);

  auto HalfEdgeIndex = [&mesh](HalfEdgeIter iter) -> int {
    return static_cast<int>(iter - mesh.halfEdges.begin());
  };

  auto MarkFlip = [&halfEdgeHasFlip, &HalfEdgeIndex](HalfEdgeIter iter, bool value) {
    const int index = HalfEdgeIndex(iter);
    if (index < 0)
    {
      return;
    }
    if (index >= static_cast<int>(halfEdgeHasFlip.size()))
    {
      halfEdgeHasFlip.resize(static_cast<std::size_t>(index) + 1, 0);
    }
    halfEdgeHasFlip[static_cast<std::size_t>(index)] = value ? 1 : 0;
  };

  auto HasFlip = [&halfEdgeHasFlip, &HalfEdgeIndex](HalfEdgeIter iter) -> bool {
    const int index = HalfEdgeIndex(iter);
    return index >= 0 && index < static_cast<int>(halfEdgeHasFlip.size()) &&
      halfEdgeHasFlip[static_cast<std::size_t>(index)] != 0;
  };

  preallocateMeshElements(data, mesh);

  for (std::size_t i = 0; i < data.positions.size(); i++)
  {
    VertexIter vertex = mesh.vertices.insert(mesh.vertices.end(), Vertex());
    vertex->position = data.positions[i];
    vertex->he = isolated.begin();
    indexToVertex[static_cast<int>(i)] = vertex;
  }

  int faceIndex = 0;
  bool hasValidFace = false;
  for (const auto& face : data.indices)
  {
    int n = static_cast<int>(face.size());

    if (n < 3)
    {
      std::cerr << "Error: face " << faceIndex << " is degenerate" << std::endl;
      faceIndex++;
      continue;
    }

    FaceIter newFace = mesh.faces.insert(mesh.faces.end(), Face());
    hasValidFace = true;

    std::vector<HalfEdgeIter> halfEdges(n);
    for (int i = 0; i < n; i++)
    {
      halfEdges[i] = mesh.halfEdges.insert(mesh.halfEdges.end(), HalfEdge());
    }

    for (int i = 0; i < n; i++)
    {
      int a = face[i].position;
      int b = face[(i + 1) % n].position;

      halfEdges[i]->next = halfEdges[(i + 1) % n];
      halfEdges[i]->vertex = indexToVertex[a];
      const int currentIndex = HalfEdgeIndex(halfEdges[i]);

      halfEdges[i]->onBoundary = false;
      MarkFlip(halfEdges[i], false);

      indexToVertex[a]->he = halfEdges[i];

      halfEdges[i]->face = newFace;
      newFace->he = halfEdges[i];

      EdgeKey edgeKey = MakeEdgeKey(a, b);
      auto existing = existingHalfEdges.find(edgeKey);
      if (existing != existingHalfEdges.end())
      {
        HalfEdgeIter opposite = mesh.halfEdges.begin() + existing->second;
        halfEdges[i]->flip = opposite;
        opposite->flip = halfEdges[i];
        halfEdges[i]->edge = opposite->edge;
        MarkFlip(halfEdges[i], true);
        MarkFlip(opposite, true);
      }
      else
      {
        halfEdges[i]->edge = mesh.edges.insert(mesh.edges.end(), Edge());
        halfEdges[i]->edge->he = halfEdges[i];
      }

      existingHalfEdges[edgeKey] = currentIndex;

      int& edgeUseCount = edgeCount[edgeKey];
      edgeUseCount++;
      if (edgeUseCount > 2)
      {
        std::cerr << "Error: edge " << edgeKey.A << ", " << edgeKey.B
                  << " is non manifold" << std::endl;
        return false;
      }
    }

    faceIndex++;
  }

  if (!hasValidFace)
  {
    return false;
  }

  for (HalfEdgeIter currHe = mesh.halfEdges.begin(); currHe != mesh.halfEdges.end(); currHe++)
  {
    if (!HasFlip(currHe))
    {
      FaceIter newFace = mesh.faces.insert(mesh.faces.end(), Face());

      std::vector<HalfEdgeIter> boundaryCycle;
      HalfEdgeIter he = currHe;
      do
      {
        HalfEdgeIter newHe = mesh.halfEdges.insert(mesh.halfEdges.end(), HalfEdge());
        newHe->onBoundary = true;

        he->flip = newHe;

        HalfEdgeIter nextHe = he->next;
        while (HasFlip(nextHe))
        {
          HalfEdgeIter flipped = nextHe->flip;
          if (flipped->onBoundary)
          {
            nextHe = currHe;
            break;
          }

          nextHe = flipped->next;
        }

        newHe->flip = he;
        newHe->vertex = nextHe->vertex;
        newHe->edge = he->edge;
        newHe->face = newFace;

        newFace->he = newHe;

        MarkFlip(newHe, true);
        MarkFlip(he, true);

        boundaryCycle.push_back(newHe);

        he = nextHe;

      } while (he != currHe);

      int n = static_cast<int>(boundaryCycle.size());
      for (int i = 0; i < n; i++)
      {
        boundaryCycle[i]->next = boundaryCycle[(i + n - 1) % n];
        MarkFlip(boundaryCycle[i], true);
        MarkFlip(boundaryCycle[i]->flip, true);
      }
      mesh.boundaries.insert(mesh.boundaries.end(), boundaryCycle[0]);
    }
  }

  indexElements(mesh);
  checkIsolatedVertices(mesh);
  checkNonManifoldVertices(mesh);

  return true;
}

bool MeshIO::read(std::ifstream& in, Mesh& mesh)
{
  MeshData data;

  std::string line;
  while (getline(in, line))
  {
    std::stringstream ss(line);
    std::string token;

    ss >> token;

    if (token == "v")
    {
      double x, y, z;
      ss >> x >> y >> z;

      data.positions.emplace_back(x, y, z);
    }
    else if (token == "vt")
    {
      double u, v;
      ss >> u >> v;

      data.uvs.emplace_back(u, v, 0.0);
    }
    else if (token == "vn")
    {
      double x, y, z;
      ss >> x >> y >> z;

      data.normals.emplace_back(x, y, z);
    }
    else if (token == "f")
    {
      std::vector<Index> faceIndices;

      while (ss >> token)
      {
  Index index = parseFaceIndex(token);
        if (index.position < 0)
        {
          getline(in, line);
          std::size_t i = line.find_first_not_of("\t\n\v\f\r ");
          index = parseFaceIndex(line.substr(i));
        }

        faceIndices.push_back(index);
      }

  data.indices.push_back(faceIndices);
    }
  }

  return buildMesh(data, mesh);
}

void MeshIO::write(std::ofstream& out, const Mesh& mesh)
{
  std::unordered_map<std::string, int> vertexMap;
  std::unordered_map<std::string, int> uvMap;
  std::unordered_map<std::string, int> normalMap;

  int index = 1;
  for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++)
  {
    out << "v " << v->position.x() << " " << v->position.y() << " " << v->position.z()
        << std::endl;

    vertexMap[stringRep(v->position)] = index;
    index++;
  }

  index = 0;
  for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++)
  {
    HalfEdgeIter he = mesh.faces[index].he;

    if (he->onBoundary)
    {
      index++;
      continue;
    }

    out << "f";
    HalfEdgeIter current = he;
    do
    {
  out << " " << vertexMap[stringRep(current->vertex->position)];
      current = current->next;

    } while (current != he);
    out << std::endl;

    index++;
  }
}

bool MeshIO::BuildFromPolyData(vtkPolyData* poly, Mesh& mesh)
{
  if (!poly)
  {
    return false;
  }

  vtkPoints* points = poly->GetPoints();
  if (!points)
  {
    return false;
  }

  MeshData data;
  vtkIdType numPoints = points->GetNumberOfPoints();
  data.positions.resize(static_cast<std::size_t>(numPoints));
  for (vtkIdType i = 0; i < numPoints; ++i)
  {
    double p[3];
    points->GetPoint(i, p);
    data.positions[static_cast<std::size_t>(i)] = Eigen::Vector3d(p[0], p[1], p[2]);
  }

  vtkCellArray* polys = poly->GetPolys();
  vtkCellArray* strips = poly->GetStrips();

  auto appendFromCellArray = [&data](vtkCellArray* cellArray) {
    if (!cellArray)
    {
      return;
    }

    cellArray->InitTraversal();
    vtkIdType npts = 0;
    const vtkIdType* pts = nullptr;
    while (cellArray->GetNextCell(npts, pts))
    {
      if (npts < 3)
      {
        continue;
      }

  std::vector<Index> faceIndices;
      faceIndices.reserve(static_cast<std::size_t>(npts));
      for (vtkIdType j = 0; j < npts; ++j)
      {
        faceIndices.emplace_back(static_cast<int>(pts[j]), -1, -1);
      }
      data.indices.push_back(std::move(faceIndices));
    }
  };

  appendFromCellArray(polys);
  appendFromCellArray(strips);

  if (data.indices.empty())
  {
    return false;
  }

  return buildMesh(data, mesh);
}

void MeshIO::CopyToPolyData(const Mesh& mesh, vtkPolyData* poly)
{
  if (!poly)
  {
    return;
  }

  vtkNew<vtkPoints> points;
  points->SetNumberOfPoints(static_cast<vtkIdType>(mesh.vertices.size()));
  for (const Vertex& v : mesh.vertices)
  {
    points->SetPoint(v.index, v.position.x(), v.position.y(), v.position.z());
  }

  vtkNew<vtkCellArray> polys;
  vtkIdType tri[3];
  for (const Face& face : mesh.faces)
  {
    if (face.remove || face.isBoundary())
    {
      continue;
    }

    HalfEdgeIter he = face.he;
    tri[0] = he->vertex->index;
    tri[1] = he->next->vertex->index;
    tri[2] = he->next->next->vertex->index;
    polys->InsertNextCell(3, tri);
  }

  poly->SetPoints(points);
  poly->SetPolys(polys);
  poly->SetVerts(nullptr);
  poly->SetLines(nullptr);
  poly->SetStrips(nullptr);
}

} // namespace vtkBotschKobbeltRemeshing
