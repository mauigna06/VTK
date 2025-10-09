#include "vtkIsotropicRemeshing.h"

#include "vtkCellArray.h"
#include "vtkIdList.h"
#include "vtkMath.h"
#include "vtkMutableDirectedGraph.h"
#include "vtkNew.h"
#include "vtkPolyData.h"
#include "vtkPoints.h"
#include "vtkTriangle.h"
#include "vtkCellData.h"
#include "vtkPointData.h"
#include "vtkSmartPointer.h"
#include "vtkObjectFactory.h"
#include "vtkDoubleArray.h"
#include "vtkCellArray.h"

// Needed for RequestData: information objects and data object keys
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkDataObject.h"

// STL
#include <vector>
#include <array>
#include <map>

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
// compute triangle area and centroid
inline double TriangleArea(const double a[3], const double b[3], const double c[3])
{
  double v1[3], v2[3], cross[3];
  vtkMath::Subtract(b, a, v1);
  vtkMath::Subtract(c, a, v2);
  vtkMath::Cross(v1, v2, cross);
  return 0.5 * vtkMath::Norm(cross);
}

inline void AddWeighted(double out[3], const double in[3], double w)
{
  out[0] += in[0] * w;
  out[1] += in[1] * w;
  out[2] += in[2] * w;
}
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

  input->BuildLinks();
  vtkNew<vtkPoints> points;
  points->DeepCopy(input->GetPoints());

  vtkIdType numPts = points->GetNumberOfPoints();
  vtkIdType numCells = input->GetNumberOfCells();

  // Precompute triangle centroids and areas
  std::vector<std::array<double,3>> triCentroid(numCells);
  std::vector<double> triArea(numCells, 0.0);

  vtkNew<vtkIdList> cellPts;
  for (vtkIdType cid = 0; cid < numCells; ++cid)
  {
    input->GetCellPoints(cid, cellPts);
    if (cellPts->GetNumberOfIds() != 3)
    {
      // skip non-triangles
      triCentroid[cid] = {0.0,0.0,0.0};
      triArea[cid] = 0.0;
      continue;
    }
    double p0[3], p1[3], p2[3];
    points->GetPoint(cellPts->GetId(0), p0);
    points->GetPoint(cellPts->GetId(1), p1);
    points->GetPoint(cellPts->GetId(2), p2);
    triArea[cid] = TriangleArea(p0,p1,p2);
    triCentroid[cid][0] = (p0[0] + p1[0] + p2[0]) / 3.0;
    triCentroid[cid][1] = (p0[1] + p1[1] + p2[1]) / 3.0;
    triCentroid[cid][2] = (p0[2] + p1[2] + p2[2]) / 3.0;
  }

  // Identify boundary points if needed
  std::vector<char> isBoundary(numPts, 0);
  if (this->PreserveTopology)
  {
    for (vtkIdType pid = 0; pid < numPts; ++pid)
    {
      // Using link-based edge counting: a boundary edge belongs to only one polygon.
      // We'll look at incident cells and edges by inspecting neighbor cells around point.
      vtkNew<vtkIdList> cellIds;
      input->GetPointCells(pid, cellIds);
      std::map<std::pair<vtkIdType,vtkIdType>, int> edgeCount;
      for (vtkIdType ci = 0; ci < cellIds->GetNumberOfIds(); ++ci)
      {
        vtkIdType cid = cellIds->GetId(ci);
        input->GetCellPoints(cid, cellPts);
        for (vtkIdType k = 0; k < cellPts->GetNumberOfIds(); ++k)
        {
          vtkIdType a = cellPts->GetId(k);
          vtkIdType b = cellPts->GetId((k+1)%cellPts->GetNumberOfIds());
          if (a > b) std::swap(a,b);
          edgeCount[{a,b}]++;
        }
      }
      for (auto &kv : edgeCount)
      {
        if (kv.second == 1)
        {
          // if any incident edge is boundary, mark point boundary
          isBoundary[pid] = 1;
          break;
        }
      }
    }
  }

  // One-ring connectivity: for each point, collect incident triangle ids
  std::vector<std::vector<vtkIdType>> incidentTris(numPts);
  for (vtkIdType cid = 0; cid < numCells; ++cid)
  {
    input->GetCellPoints(cid, cellPts);
    for (vtkIdType k = 0; k < cellPts->GetNumberOfIds(); ++k)
    {
      vtkIdType pid = cellPts->GetId(k);
      incidentTris[pid].push_back(cid);
    }
  }

  // Iterative Lloyd-like relaxation: move vertices to area-weighted centroid of triangle centroids
  std::vector<double> newPos(3*numPts);
  for (int it = 0; it < this->NumberOfIterations; ++it)
  {
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

      double accum[3] = {0.0,0.0,0.0};
      double wsum = 0.0;
      for (vtkIdType cid : incidentTris[pid])
      {
        double area = triArea[cid];
        if (area <= 0.0) continue;
        AddWeighted(accum, triCentroid[cid].data(), area);
        wsum += area;
      }
      if (wsum > 0.0)
      {
        double centroid[3] = {accum[0]/wsum, accum[1]/wsum, accum[2]/wsum};
        // move towards centroid by relaxation factor
        newPos[3*pid+0] = p[0] + this->RelaxationFactor * (centroid[0] - p[0]);
        newPos[3*pid+1] = p[1] + this->RelaxationFactor * (centroid[1] - p[1]);
        newPos[3*pid+2] = p[2] + this->RelaxationFactor * (centroid[2] - p[2]);
      }
      else
      {
        newPos[3*pid+0] = p[0];
        newPos[3*pid+1] = p[1];
        newPos[3*pid+2] = p[2];
      }
    }

    // Apply new positions
    for (vtkIdType pid = 0; pid < numPts; ++pid)
    {
      points->SetPoint(pid, &newPos[3*pid]);
    }

    // Recompute triangle centroids (areas unchanged if topology unchanged)
    for (vtkIdType cid = 0; cid < numCells; ++cid)
    {
      input->GetCellPoints(cid, cellPts);
      if (cellPts->GetNumberOfIds() != 3) continue;
      double p0[3], p1[3], p2[3];
      points->GetPoint(cellPts->GetId(0), p0);
      points->GetPoint(cellPts->GetId(1), p1);
      points->GetPoint(cellPts->GetId(2), p2);
      triCentroid[cid][0] = (p0[0] + p1[0] + p2[0]) / 3.0;
      triCentroid[cid][1] = (p0[1] + p1[1] + p2[1]) / 3.0;
      triCentroid[cid][2] = (p0[2] + p1[2] + p2[2]) / 3.0;
      // area may be recomputed if desired; we keep it from original mesh to stabilize weights
    }
  }

  // Output: copy input topology and use new points
  output->ShallowCopy(input);
  output->SetPoints(points);

  return 1;
}
