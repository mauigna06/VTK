#include "vtkIsotropicRemeshing.h"

#include "vtkIdList.h"
#include "vtkMath.h"
#include "vtkNew.h"
#include "vtkObjectFactory.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"

// Needed for RequestData: information objects and data object keys
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkDataObject.h"

// STL
#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <unordered_map>
#include <vector>

vtkStandardNewMacro(vtkIsotropicRemeshing);

vtkIsotropicRemeshing::vtkIsotropicRemeshing()
{
  this->NumberOfIterations = 10;
  this->RelaxationFactor = 1.0;
  this->TargetEdgeLength = -1.0;
  this->PreserveTopology = true;
}

vtkIsotropicRemeshing::~vtkIsotropicRemeshing() = default;

void vtkIsotropicRemeshing::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
  os << indent << "NumberOfIterations: " << this->NumberOfIterations << "\n";
  os << indent << "RelaxationFactor: " << this->RelaxationFactor << "\n";
  os << indent << "TargetEdgeLength: " << this->TargetEdgeLength << "\n";
  os << indent << "PreserveTopology: " << this->PreserveTopology << "\n";
}

namespace
{
struct EdgeKey
{
  vtkIdType A;
  vtkIdType B;

  bool operator==(const EdgeKey& other) const noexcept
  {
    return this->A == other.A && this->B == other.B;
  }
};

struct EdgeKeyHash
{
  std::size_t operator()(const EdgeKey& e) const noexcept
  {
    std::size_t h1 = std::hash<vtkIdType>{}(e.A);
    std::size_t h2 = std::hash<vtkIdType>{}(e.B);
    return h1 ^ (h2 + 0x9e3779b97f4a7c15ULL + (h1 << 6) + (h1 >> 2));
  }
};
}

int vtkIsotropicRemeshing::RequestData(vtkInformation* vtkNotUsed(request),
                                        vtkInformationVector** inputVector,
                                        vtkInformationVector* outputVector)
{
  vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
  vtkPolyData* input = vtkPolyData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkInformation* outInfo = outputVector->GetInformationObject(0);
  vtkPolyData* output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

  if (!input || !output)
  {
    vtkErrorMacro("Invalid input or output");
    return 0;
  }

  vtkNew<vtkPoints> points;
  points->DeepCopy(input->GetPoints());

  vtkIdType numPts = points->GetNumberOfPoints();
  vtkIdType numCells = input->GetNumberOfCells();

  if (numPts == 0 || numCells == 0)
  {
    output->ShallowCopy(input);
    return 1;
  }
  vtkNew<vtkIdList> cellPts;

  // Build adjacency and edge statistics
  std::vector<std::vector<vtkIdType>> incidentTris(numPts);
  std::vector<std::vector<vtkIdType>> vertexNeighbors(numPts);
  std::unordered_map<EdgeKey, int, EdgeKeyHash> edgeUse;
  edgeUse.reserve(static_cast<size_t>(numCells) * 3);

  for (vtkIdType cid = 0; cid < numCells; ++cid)
  {
    input->GetCellPoints(cid, cellPts);
    if (cellPts->GetNumberOfIds() != 3)
    {
      continue;
    }

    vtkIdType id0 = cellPts->GetId(0);
    vtkIdType id1 = cellPts->GetId(1);
    vtkIdType id2 = cellPts->GetId(2);

    incidentTris[id0].push_back(cid);
    incidentTris[id1].push_back(cid);
    incidentTris[id2].push_back(cid);

    vertexNeighbors[id0].push_back(id1);
    vertexNeighbors[id0].push_back(id2);
    vertexNeighbors[id1].push_back(id0);
    vertexNeighbors[id1].push_back(id2);
    vertexNeighbors[id2].push_back(id0);
    vertexNeighbors[id2].push_back(id1);

    EdgeKey e0{ std::min(id0, id1), std::max(id0, id1) };
    EdgeKey e1{ std::min(id1, id2), std::max(id1, id2) };
    EdgeKey e2{ std::min(id2, id0), std::max(id2, id0) };
    edgeUse[e0]++;
    edgeUse[e1]++;
    edgeUse[e2]++;
  }

  for (auto& nbs : vertexNeighbors)
  {
    std::sort(nbs.begin(), nbs.end());
    nbs.erase(std::unique(nbs.begin(), nbs.end()), nbs.end());
  }

  std::vector<EdgeKey> edges;
  edges.reserve(edgeUse.size());
  for (const auto& kv : edgeUse)
  {
    edges.push_back(kv.first);
  }

  // Identify boundary vertices if topology preservation is requested
  std::vector<char> isBoundary(numPts, 0);
  if (this->PreserveTopology)
  {
    for (const auto& kv : edgeUse)
    {
      if (kv.second == 1)
      {
        isBoundary[kv.first.A] = 1;
        isBoundary[kv.first.B] = 1;
      }
    }
  }

  // Estimate target edge length if not provided
  double targetEdgeLength = this->TargetEdgeLength;
  if (targetEdgeLength <= 0.0 && !edgeUse.empty())
  {
    double totalLength = 0.0;
    vtkIdType edgeCount = 0;
    double pa[3], pb[3];
    for (const auto& kv : edgeUse)
    {
      points->GetPoint(kv.first.A, pa);
      points->GetPoint(kv.first.B, pb);
      double dist = std::sqrt(vtkMath::Distance2BetweenPoints(pa, pb));
      if (dist > 0.0)
      {
        totalLength += dist;
        ++edgeCount;
      }
    }
    if (edgeCount > 0)
    {
      targetEdgeLength = totalLength / static_cast<double>(edgeCount);
    }
  }

  const double localEdgeBalancingFactor = 0.5;
  const double globalEdgeBalancingGain = 0.8;
  const double minEdgeLength = 1e-12;

  // Iterative Lloyd-like relaxation: move vertices to area-weighted centroid of triangle centroids
  std::vector<double> newPos(3*numPts);
  std::vector<double> vertexNormals(3 * numPts, 0.0);
  std::vector<double> edgeAccum(3 * numPts, 0.0);
  std::vector<int> edgeCounts(numPts, 0);
  for (int it = 0; it < this->NumberOfIterations; ++it)
  {
    std::fill(vertexNormals.begin(), vertexNormals.end(), 0.0);
    std::fill(edgeAccum.begin(), edgeAccum.end(), 0.0);
    std::fill(edgeCounts.begin(), edgeCounts.end(), 0);

    for (vtkIdType pid = 0; pid < numPts; ++pid)
    {
      double p[3];
      points->GetPoint(pid, p);
      if (this->PreserveTopology && isBoundary[pid])
      {
        // keep boundary points fixed
        newPos[3*pid+0] = p[0];
        newPos[3*pid+1] = p[1];
        newPos[3*pid+2] = p[2];
        continue;
      }

      double centroidAccum[3] = { 0.0, 0.0, 0.0 };
      double areaSum = 0.0;
      double normalAccum[3] = { 0.0, 0.0, 0.0 };

      for (vtkIdType cid : incidentTris[pid])
      {
        input->GetCellPoints(cid, cellPts);
        if (cellPts->GetNumberOfIds() != 3)
        {
          continue;
        }

        double p0[3], p1[3], p2[3];
        points->GetPoint(cellPts->GetId(0), p0);
        points->GetPoint(cellPts->GetId(1), p1);
        points->GetPoint(cellPts->GetId(2), p2);

        double v1[3], v2[3], triNormal[3];
        vtkMath::Subtract(p1, p0, v1);
        vtkMath::Subtract(p2, p0, v2);
        vtkMath::Cross(v1, v2, triNormal);
        double area = 0.5 * vtkMath::Norm(triNormal);
        if (area <= 0.0)
        {
          continue;
        }

        double centroidTri[3] = { (p0[0] + p1[0] + p2[0]) / 3.0,
          (p0[1] + p1[1] + p2[1]) / 3.0, (p0[2] + p1[2] + p2[2]) / 3.0 };
        centroidAccum[0] += centroidTri[0] * area;
        centroidAccum[1] += centroidTri[1] * area;
        centroidAccum[2] += centroidTri[2] * area;
        areaSum += area;

        normalAccum[0] += triNormal[0];
        normalAccum[1] += triNormal[1];
        normalAccum[2] += triNormal[2];
      }

      double centroidMove[3] = { 0.0, 0.0, 0.0 };
      if (areaSum > 0.0)
      {
        double inv = 1.0 / areaSum;
        double centroid[3] = { centroidAccum[0] * inv, centroidAccum[1] * inv,
          centroidAccum[2] * inv };
        centroidMove[0] = centroid[0] - p[0];
        centroidMove[1] = centroid[1] - p[1];
        centroidMove[2] = centroid[2] - p[2];
      }

      double normalLen = vtkMath::Norm(normalAccum);
      double normal[3] = { 0.0, 0.0, 0.0 };
      if (normalLen > 0.0)
      {
        normal[0] = normalAccum[0] / normalLen;
        normal[1] = normalAccum[1] / normalLen;
        normal[2] = normalAccum[2] / normalLen;

        double dot = vtkMath::Dot(centroidMove, normal);
        centroidMove[0] -= dot * normal[0];
        centroidMove[1] -= dot * normal[1];
        centroidMove[2] -= dot * normal[2];
      }

      vertexNormals[3 * pid + 0] = normal[0];
      vertexNormals[3 * pid + 1] = normal[1];
      vertexNormals[3 * pid + 2] = normal[2];

      double edgeAdjust[3] = { 0.0, 0.0, 0.0 };
      vtkIdType neighborCount = 0;
      if (targetEdgeLength > 0.0)
      {
        for (vtkIdType nb : vertexNeighbors[pid])
        {
          double q[3];
          points->GetPoint(nb, q);
          double diff[3];
          vtkMath::Subtract(p, q, diff);
          double len = vtkMath::Norm(diff);
          if (len <= minEdgeLength)
          {
            continue;
          }
          double diffScale = (len - targetEdgeLength) / len;
          edgeAdjust[0] += diff[0] * diffScale;
          edgeAdjust[1] += diff[1] * diffScale;
          edgeAdjust[2] += diff[2] * diffScale;
          ++neighborCount;
        }
        if (neighborCount > 0)
        {
          double inv = 1.0 / static_cast<double>(neighborCount);
          edgeAdjust[0] *= inv;
          edgeAdjust[1] *= inv;
          edgeAdjust[2] *= inv;
          if (normalLen > 0.0)
          {
            double dot = vtkMath::Dot(edgeAdjust, normal);
            edgeAdjust[0] -= dot * normal[0];
            edgeAdjust[1] -= dot * normal[1];
            edgeAdjust[2] -= dot * normal[2];
          }
        }
      }

      double proposed[3] = { p[0], p[1], p[2] };
      proposed[0] += this->RelaxationFactor * centroidMove[0];
      proposed[1] += this->RelaxationFactor * centroidMove[1];
      proposed[2] += this->RelaxationFactor * centroidMove[2];

      if (targetEdgeLength > 0.0 && neighborCount > 0)
      {
        double scale = this->RelaxationFactor * localEdgeBalancingFactor;
        proposed[0] -= scale * edgeAdjust[0];
        proposed[1] -= scale * edgeAdjust[1];
        proposed[2] -= scale * edgeAdjust[2];
      }

      double move[3] = { proposed[0] - p[0], proposed[1] - p[1], proposed[2] - p[2] };
      double moveLen = vtkMath::Norm(move);
      if (targetEdgeLength > 0.0)
      {
        double maxStep = 0.5 * targetEdgeLength;
        if (moveLen > maxStep && maxStep > 0.0)
        {
          double s = maxStep / moveLen;
          move[0] *= s;
          move[1] *= s;
          move[2] *= s;
          proposed[0] = p[0] + move[0];
          proposed[1] = p[1] + move[1];
          proposed[2] = p[2] + move[2];
        }
      }

      newPos[3 * pid + 0] = proposed[0];
      newPos[3 * pid + 1] = proposed[1];
      newPos[3 * pid + 2] = proposed[2];
    }

    double lengthSum = 0.0;
    vtkIdType lengthCount = 0;

    for (const EdgeKey& edge : edges)
    {
      vtkIdType a = edge.A;
      vtkIdType b = edge.B;
      if (a < 0 || b < 0 || a >= numPts || b >= numPts)
      {
        continue;
      }

      const double* aPos = &newPos[3 * a];
      const double* bPos = &newPos[3 * b];
      double diff[3] = { bPos[0] - aPos[0], bPos[1] - aPos[1], bPos[2] - aPos[2] };
      double len = vtkMath::Norm(diff);
      if (len <= minEdgeLength)
      {
        continue;
      }

      lengthSum += len;
      ++lengthCount;

      double unit[3] = { diff[0] / len, diff[1] / len, diff[2] / len };
      double delta = 0.5 * (len - targetEdgeLength);

      if (!(this->PreserveTopology && isBoundary[a]))
      {
        double corr[3] = { unit[0] * delta, unit[1] * delta, unit[2] * delta };
        double normal[3] = { vertexNormals[3 * a + 0], vertexNormals[3 * a + 1],
          vertexNormals[3 * a + 2] };
        double nLen = vtkMath::Norm(normal);
        if (nLen > 0.0)
        {
          double dot = vtkMath::Dot(corr, normal);
          corr[0] -= dot * normal[0];
          corr[1] -= dot * normal[1];
          corr[2] -= dot * normal[2];
        }
        edgeAccum[3 * a + 0] += corr[0];
        edgeAccum[3 * a + 1] += corr[1];
        edgeAccum[3 * a + 2] += corr[2];
        ++edgeCounts[a];
      }

      if (!(this->PreserveTopology && isBoundary[b]))
      {
        double corr[3] = { -unit[0] * delta, -unit[1] * delta, -unit[2] * delta };
        double normal[3] = { vertexNormals[3 * b + 0], vertexNormals[3 * b + 1],
          vertexNormals[3 * b + 2] };
        double nLen = vtkMath::Norm(normal);
        if (nLen > 0.0)
        {
          double dot = vtkMath::Dot(corr, normal);
          corr[0] -= dot * normal[0];
          corr[1] -= dot * normal[1];
          corr[2] -= dot * normal[2];
        }
        edgeAccum[3 * b + 0] += corr[0];
        edgeAccum[3 * b + 1] += corr[1];
        edgeAccum[3 * b + 2] += corr[2];
        ++edgeCounts[b];
      }
    }

    if (this->TargetEdgeLength <= 0.0 && lengthCount > 0)
    {
      targetEdgeLength = lengthSum / static_cast<double>(lengthCount);
    }

    for (vtkIdType pid = 0; pid < numPts; ++pid)
    {
      if (this->PreserveTopology && isBoundary[pid])
      {
        continue;
      }

      if (edgeCounts[pid] > 0)
      {
        double corr[3] = { edgeAccum[3 * pid + 0] / static_cast<double>(edgeCounts[pid]),
          edgeAccum[3 * pid + 1] / static_cast<double>(edgeCounts[pid]),
          edgeAccum[3 * pid + 2] / static_cast<double>(edgeCounts[pid]) };

        double normal[3] = { vertexNormals[3 * pid + 0], vertexNormals[3 * pid + 1],
          vertexNormals[3 * pid + 2] };
        double nLen = vtkMath::Norm(normal);
        if (nLen > 0.0)
        {
          double dot = vtkMath::Dot(corr, normal);
          corr[0] -= dot * normal[0];
          corr[1] -= dot * normal[1];
          corr[2] -= dot * normal[2];
        }

        double corrLen = vtkMath::Norm(corr);
        if (targetEdgeLength > 0.0)
        {
          double clamp = 0.25 * targetEdgeLength;
          if (clamp > 0.0 && corrLen > clamp)
          {
            double s = clamp / corrLen;
            corr[0] *= s;
            corr[1] *= s;
            corr[2] *= s;
          }
        }

        newPos[3 * pid + 0] += globalEdgeBalancingGain * this->RelaxationFactor * corr[0];
        newPos[3 * pid + 1] += globalEdgeBalancingGain * this->RelaxationFactor * corr[1];
        newPos[3 * pid + 2] += globalEdgeBalancingGain * this->RelaxationFactor * corr[2];
      }
    }

    // Apply new positions
    for (vtkIdType pid = 0; pid < numPts; ++pid)
    {
      points->SetPoint(pid, &newPos[3*pid]);
    }

    points->Modified();
  }

  // Output: copy input topology and use new points
  output->ShallowCopy(input);
  output->SetPoints(points);

  return 1;
}
