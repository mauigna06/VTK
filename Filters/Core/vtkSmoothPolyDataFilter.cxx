// SPDX-FileCopyrightText: Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
// SPDX-License-Identifier: BSD-3-Clause
#include "vtkSmoothPolyDataFilter.h"

#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkCellLocator.h"
#include "vtkFloatArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkMath.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkPolygon.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkTriangleFilter.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <memory>
#include <unordered_set>
#include <vector>

VTK_ABI_NAMESPACE_BEGIN
vtkStandardNewMacro(vtkSmoothPolyDataFilter);

// The following code defines a helper class for performing mesh smoothing
// across the surface of another mesh.
typedef struct
{
  vtkIdType cellId; // cell
  int subId;        // cell sub id
  double p[3];      // parametric coords in cell
} vtkSmoothPoint;

class vtkSmoothPoints
{ //;prevent man page generation
public:
  vtkSmoothPoints();
  ~vtkSmoothPoints() { delete[] this->Array; }
  vtkIdType GetNumberOfPoints() { return this->MaxId + 1; }
  vtkSmoothPoint* GetSmoothPoint(vtkIdType i) { return this->Array + i; }
  vtkSmoothPoint* InsertSmoothPoint(vtkIdType ptId)
  {
    if (ptId >= this->Size)
    {
      this->Resize(ptId + 1);
    }
    if (ptId > this->MaxId)
    {
      this->MaxId = ptId;
    }
    return this->Array + ptId;
  }
  vtkSmoothPoint* Resize(vtkIdType sz); // reallocates data
  void Reset() { this->MaxId = -1; }

  vtkSmoothPoint* Array; // pointer to data
  vtkIdType MaxId;       // maximum index inserted thus far
  vtkIdType Size;        // allocated size of data
  vtkIdType Extend;      // grow array by this amount
};

//------------------------------------------------------------------------------
vtkSmoothPoints::vtkSmoothPoints()
{
  this->MaxId = -1;
  this->Array = new vtkSmoothPoint[1000];
  this->Size = 1000;
  this->Extend = 5000;
}

//------------------------------------------------------------------------------
vtkSmoothPoint* vtkSmoothPoints::Resize(vtkIdType sz)
{
  vtkSmoothPoint* newArray;
  vtkIdType newSize;

  if (sz >= this->Size)
  {
    newSize = this->Size + this->Extend * (((sz - this->Size) / this->Extend) + 1);
  }
  else
  {
    newSize = sz;
  }

  newArray = new vtkSmoothPoint[newSize];

  memcpy(newArray, this->Array, (sz < this->Size ? sz : this->Size) * sizeof(vtkSmoothPoint));

  this->Size = newSize;
  delete[] this->Array;
  this->Array = newArray;

  return this->Array;
}

// The following code defines methods for the vtkSmoothPolyDataFilter class
//

// Construct object with number of iterations 20; relaxation factor .01;
// feature edge smoothing turned off; feature
// angle 45 degrees; edge angle 15 degrees; and boundary smoothing turned
// on. Error scalars and vectors are not generated (by default). The
// convergence criterion is 0.0 of the bounding box diagonal.
//------------------------------------------------------------------------------
vtkSmoothPolyDataFilter::vtkSmoothPolyDataFilter()
{
  this->Convergence = 0.0; // goes to number of specified iterations
  this->NumberOfIterations = 20;

  this->RelaxationFactor = .01;

  this->FeatureAngle = 45.0;
  this->EdgeAngle = 15.0;
  this->FeatureEdgeSmoothing = 0;
  this->BoundarySmoothing = 1;

  this->GenerateErrorScalars = 0;
  this->GenerateErrorVectors = 0;

  this->OutputPointsPrecision = vtkAlgorithm::DEFAULT_PRECISION;
  this->SmoothingMode = VTK_SMOOTH_MODE_LAPLACIAN;

  this->SmoothPoints = nullptr;

  // optional second input
  this->SetNumberOfInputPorts(2);
}

//------------------------------------------------------------------------------
vtkSmoothPolyDataFilter::~vtkSmoothPolyDataFilter() = default;

//------------------------------------------------------------------------------
void vtkSmoothPolyDataFilter::SetSourceData(vtkPolyData* source)
{
  this->SetInputData(1, source);
}

//------------------------------------------------------------------------------
vtkPolyData* vtkSmoothPolyDataFilter::GetSource()
{
  if (this->GetNumberOfInputConnections(1) < 1)
  {
    return nullptr;
  }
  return vtkPolyData::SafeDownCast(this->GetExecutive()->GetInputData(1, 0));
}

//------------------------------------------------------------------------------
void vtkSmoothPolyDataFilter::SetUniformTriangleSmoothing(vtkTypeBool value)
{
  this->SetSmoothingMode(
    value ? VTK_SMOOTH_MODE_UNIFORM_TRIANGLE : VTK_SMOOTH_MODE_LAPLACIAN);
}

//------------------------------------------------------------------------------
vtkTypeBool vtkSmoothPolyDataFilter::GetUniformTriangleSmoothing()
{
  return (this->SmoothingMode == VTK_SMOOTH_MODE_UNIFORM_TRIANGLE) ? 1 : 0;
}

#define VTK_SIMPLE_VERTEX 0
#define VTK_FIXED_VERTEX 1
#define VTK_FEATURE_EDGE_VERTEX 2
#define VTK_BOUNDARY_EDGE_VERTEX 3

namespace
{

// Special structure for marking vertices
typedef struct _vtkMeshVertex
{
  char type;
  vtkIdList* edges; // connected edges (list of connected point ids)
  _vtkMeshVertex()
  {
    type = VTK_SIMPLE_VERTEX; // can smooth
    edges = nullptr;
  }
} vtkMeshVertex, *vtkMeshVertexPtr;

template <typename T>
struct vtkSPDF_InternalParams
{
  vtkSmoothPolyDataFilter* spdf;
  int numberOfIterations;
  vtkPoints* newPts;
  T factor;
  T conv;
  vtkIdType numPts;
  vtkMeshVertexPtr vertexPtr;
  vtkPolyData* source;
  vtkSmoothPoints* SmoothPoints;
  double* w;
  vtkCellLocator* cellLocator;
};

template <typename T>
void vtkSPDF_MovePoints(vtkSPDF_InternalParams<T>& params)
{
  int iterationNumber = 0;
  for (T maxDist = std::numeric_limits<T>::max();
       maxDist > params.conv && iterationNumber < params.numberOfIterations; ++iterationNumber)
  {
    if (iterationNumber && !(iterationNumber % 5))
    {
      params.spdf->UpdateProgress(0.5 + 0.5 * iterationNumber / params.numberOfIterations);
      if (params.spdf->CheckAbort())
      {
        break;
      }
    }

    maxDist = 0.0;
    T* newPtsCoords = static_cast<T*>(params.newPts->GetVoidPointer(0));
    T* start = newPtsCoords;
    vtkMeshVertexPtr vertsPtr = params.vertexPtr;
    vtkIdType npts, *edgeIdPtr;
    T dist, deltaX[3];
    double dist2, xNew[3], closestPt[3];

    // For each non-fixed vertex of the mesh, move the point toward the mean
    // position of its connected neighbors using the relaxation factor.
    for (vtkIdType i = 0; i < params.numPts; ++i)
    {
      if (vertsPtr->type != VTK_FIXED_VERTEX && vertsPtr->edges &&
        (npts = vertsPtr->edges->GetNumberOfIds()) > 0)
      {
        deltaX[0] = deltaX[1] = deltaX[2] = 0.0;
        edgeIdPtr = vertsPtr->edges->GetPointer(0);
        // Compute the mean (cumulated) direction vector
        for (vtkIdType j = 0; j < npts; ++j)
        {
          for (unsigned short k = 0; k < 3; ++k)
          {
            deltaX[k] += *(start + 3 * (*edgeIdPtr) + k);
          }
          ++edgeIdPtr;
        } // for all connected points

        // Move the point
        *newPtsCoords += params.factor * (deltaX[0] / npts - (*newPtsCoords));
        xNew[0] = *newPtsCoords;
        ++newPtsCoords;
        *newPtsCoords += params.factor * (deltaX[1] / npts - (*newPtsCoords));
        xNew[1] = *newPtsCoords;
        ++newPtsCoords;
        *newPtsCoords += params.factor * (deltaX[2] / npts - (*newPtsCoords));
        xNew[2] = *newPtsCoords;
        ++newPtsCoords;

        // Constrain point to surface
        if (params.source)
        {
          vtkSmoothPoint* sPtr = params.SmoothPoints->GetSmoothPoint(i);
          vtkCell* cell = nullptr;

          if (sPtr->cellId >= 0) // in cell
          {
            cell = params.source->GetCell(sPtr->cellId);
          }

          if (!cell ||
            cell->EvaluatePosition(xNew, closestPt, sPtr->subId, sPtr->p, dist2, params.w) == 0)
          { // not in cell anymore
            params.cellLocator->FindClosestPoint(xNew, closestPt, sPtr->cellId, sPtr->subId, dist2);
          }
          for (int k = 0; k < 3; ++k)
          {
            xNew[k] = closestPt[k];
          }
          params.newPts->SetPoint(i, xNew);
        }

        if ((dist = vtkMath::Norm(deltaX)) > maxDist)
        {
          maxDist = dist;
        }
      } // if can move point
      else
      {
        newPtsCoords += 3;
      }
      ++vertsPtr;
    } // for all points
  }   // for not converged or within iteration count

  vtkDebugWithObjectMacro(params.spdf, << "Performed " << iterationNumber << " smoothing passes");
}

struct vtkSPDF_EdgeKey
{
  vtkIdType A;
  vtkIdType B;

  bool operator==(const vtkSPDF_EdgeKey& other) const noexcept
  {
    return this->A == other.A && this->B == other.B;
  }
};

struct vtkSPDF_EdgeKeyHash
{
  std::size_t operator()(const vtkSPDF_EdgeKey& e) const noexcept
  {
    std::size_t h1 = std::hash<vtkIdType>{}(e.A);
    std::size_t h2 = std::hash<vtkIdType>{}(e.B);
    return h1 ^ (h2 + 0x9e3779b97f4a7c15ULL + (h1 << 6) + (h1 >> 2));
  }
};

template <typename T>
void vtkSPDF_MovePointsUniformTriangles(vtkSPDF_InternalParams<T>& params,
  const std::vector<std::array<vtkIdType, 3>>& triangles,
  const std::vector<std::vector<vtkIdType>>& incidentTris)
{
  if (triangles.empty())
  {
    vtkSPDF_MovePoints(params);
    return;
  }

  std::unordered_set<vtkSPDF_EdgeKey, vtkSPDF_EdgeKeyHash> edgeSet;
  edgeSet.reserve(triangles.size() * 3);
  std::vector<vtkSPDF_EdgeKey> edges;
  edges.reserve(triangles.size() * 3);

  for (const auto& tri : triangles)
  {
    vtkSPDF_EdgeKey e0{ std::min(tri[0], tri[1]), std::max(tri[0], tri[1]) };
    vtkSPDF_EdgeKey e1{ std::min(tri[1], tri[2]), std::max(tri[1], tri[2]) };
    vtkSPDF_EdgeKey e2{ std::min(tri[2], tri[0]), std::max(tri[2], tri[0]) };

    if (edgeSet.insert(e0).second)
    {
      edges.push_back(e0);
    }
    if (edgeSet.insert(e1).second)
    {
      edges.push_back(e1);
    }
    if (edgeSet.insert(e2).second)
    {
      edges.push_back(e2);
    }
  }

  std::vector<double> coords(static_cast<size_t>(params.numPts) * 3);
  auto copyCurrentPoints = [&coords, &params]() {
    T* raw = static_cast<T*>(params.newPts->GetVoidPointer(0));
    const size_t total = coords.size();
    for (size_t idx = 0; idx < total; ++idx)
    {
      coords[idx] = static_cast<double>(raw[idx]);
    }
  };

  copyCurrentPoints();

  double targetEdgeLength = 0.0;
  if (!edges.empty())
  {
    double totalLength = 0.0;
    vtkIdType count = 0;
    for (const auto& edge : edges)
    {
      const double* a = &coords[3 * static_cast<size_t>(edge.A)];
      const double* b = &coords[3 * static_cast<size_t>(edge.B)];
      double diff[3] = { a[0] - b[0], a[1] - b[1], a[2] - b[2] };
      double len = vtkMath::Norm(diff);
      if (len > 0.0)
      {
        totalLength += len;
        ++count;
      }
    }
    if (count > 0)
    {
      targetEdgeLength = totalLength / static_cast<double>(count);
    }
  }

  const double convergence = static_cast<double>(params.conv);
  const double relaxation = static_cast<double>(params.factor);
  const double edgeBalancingFactor = 0.5;
  const double clampRatio = 0.5;
  const double minEdgeLength = 1e-12;

  int iterationNumber = 0;
  double maxDist = std::numeric_limits<double>::max();
  while (maxDist > convergence && iterationNumber < params.numberOfIterations)
  {
    if (iterationNumber && !(iterationNumber % 5))
    {
      params.spdf->UpdateProgress(0.5 + 0.5 * iterationNumber / params.numberOfIterations);
      if (params.spdf->CheckAbort())
      {
        break;
      }
    }

    std::vector<std::array<double, 3>> triCentroids(triangles.size());
    std::vector<double> triAreas(triangles.size());
    std::vector<std::array<double, 3>> triNormals(triangles.size());

    for (size_t tid = 0; tid < triangles.size(); ++tid)
    {
      const auto& tri = triangles[tid];
      const double* a = &coords[3 * static_cast<size_t>(tri[0])];
      const double* b = &coords[3 * static_cast<size_t>(tri[1])];
      const double* c = &coords[3 * static_cast<size_t>(tri[2])];

      double ab[3] = { b[0] - a[0], b[1] - a[1], b[2] - a[2] };
      double ac[3] = { c[0] - a[0], c[1] - a[1], c[2] - a[2] };
      double normal[3];
      vtkMath::Cross(ab, ac, normal);
      double area = 0.5 * vtkMath::Norm(normal);

      triAreas[tid] = area;
      triCentroids[tid] = { (a[0] + b[0] + c[0]) / 3.0, (a[1] + b[1] + c[1]) / 3.0,
        (a[2] + b[2] + c[2]) / 3.0 };
      triNormals[tid] = { normal[0], normal[1], normal[2] };
    }

    maxDist = 0.0;

    for (vtkIdType pid = 0; pid < params.numPts; ++pid)
    {
      vtkMeshVertexPtr vertex = params.vertexPtr + pid;
      vtkIdList* edgeList = vertex->edges;
      vtkIdType edgeCount = edgeList ? edgeList->GetNumberOfIds() : 0;

      const double* prevPt = &coords[3 * static_cast<size_t>(pid)];
      double proposed[3] = { prevPt[0], prevPt[1], prevPt[2] };

      if (vertex->type != VTK_FIXED_VERTEX && edgeCount > 0)
      {
        const std::vector<vtkIdType>& incident = incidentTris[static_cast<size_t>(pid)];
        double centroidAccum[3] = { 0.0, 0.0, 0.0 };
        double normalAccum[3] = { 0.0, 0.0, 0.0 };
        double areaSum = 0.0;

        for (vtkIdType triId : incident)
        {
          if (triId < 0 || static_cast<size_t>(triId) >= triAreas.size())
          {
            continue;
          }
          double area = triAreas[static_cast<size_t>(triId)];
          if (area <= 0.0)
          {
            continue;
          }

          const auto& centroid = triCentroids[static_cast<size_t>(triId)];
          centroidAccum[0] += centroid[0] * area;
          centroidAccum[1] += centroid[1] * area;
          centroidAccum[2] += centroid[2] * area;

          const auto& normal = triNormals[static_cast<size_t>(triId)];
          normalAccum[0] += normal[0];
          normalAccum[1] += normal[1];
          normalAccum[2] += normal[2];

          areaSum += area;
        }

        double tangentMove[3] = { 0.0, 0.0, 0.0 };
        double normalLen = vtkMath::Norm(normalAccum);

        if (areaSum > 0.0)
        {
          double invArea = 1.0 / areaSum;
          double centroid[3] = { centroidAccum[0] * invArea, centroidAccum[1] * invArea,
            centroidAccum[2] * invArea };
          double moveVec[3] = { centroid[0] - prevPt[0], centroid[1] - prevPt[1],
            centroid[2] - prevPt[2] };

          if (normalLen > 1e-12)
          {
            double nUnit[3] = { normalAccum[0] / normalLen, normalAccum[1] / normalLen,
              normalAccum[2] / normalLen };
            double proj = vtkMath::Dot(moveVec, nUnit);
            moveVec[0] -= proj * nUnit[0];
            moveVec[1] -= proj * nUnit[1];
            moveVec[2] -= proj * nUnit[2];
          }

          tangentMove[0] = moveVec[0];
          tangentMove[1] = moveVec[1];
          tangentMove[2] = moveVec[2];
        }

        std::vector<vtkIdType> neighbors;
        neighbors.reserve(static_cast<size_t>(edgeCount));
        for (vtkIdType idx = 0; idx < edgeCount; ++idx)
        {
          vtkIdType nb = edgeList->GetId(idx);
          if (nb < 0 || nb >= params.numPts)
          {
            continue;
          }
          if (std::find(neighbors.begin(), neighbors.end(), nb) == neighbors.end())
          {
            neighbors.push_back(nb);
          }
        }

        double edgeAdjust[3] = { 0.0, 0.0, 0.0 };
        if (!neighbors.empty() && targetEdgeLength > 0.0)
        {
          for (vtkIdType nb : neighbors)
          {
            const double* npt = &coords[3 * static_cast<size_t>(nb)];
            double diff[3] = { prevPt[0] - npt[0], prevPt[1] - npt[1], prevPt[2] - npt[2] };
            double len = vtkMath::Norm(diff);
            if (len <= minEdgeLength)
            {
              continue;
            }
            double scale = (len - targetEdgeLength) / len;
            edgeAdjust[0] += diff[0] * scale;
            edgeAdjust[1] += diff[1] * scale;
            edgeAdjust[2] += diff[2] * scale;
          }

          double invCount = 1.0 / static_cast<double>(neighbors.size());
          edgeAdjust[0] *= invCount;
          edgeAdjust[1] *= invCount;
          edgeAdjust[2] *= invCount;

          if (normalLen > 1e-12)
          {
            double nUnit[3] = { normalAccum[0] / normalLen, normalAccum[1] / normalLen,
              normalAccum[2] / normalLen };
            double proj = vtkMath::Dot(edgeAdjust, nUnit);
            edgeAdjust[0] -= proj * nUnit[0];
            edgeAdjust[1] -= proj * nUnit[1];
            edgeAdjust[2] -= proj * nUnit[2];
          }
        }

        proposed[0] = prevPt[0] + relaxation * tangentMove[0];
        proposed[1] = prevPt[1] + relaxation * tangentMove[1];
        proposed[2] = prevPt[2] + relaxation * tangentMove[2];

        if (!neighbors.empty() && targetEdgeLength > 0.0)
        {
          proposed[0] -= relaxation * edgeBalancingFactor * edgeAdjust[0];
          proposed[1] -= relaxation * edgeBalancingFactor * edgeAdjust[1];
          proposed[2] -= relaxation * edgeBalancingFactor * edgeAdjust[2];
        }

        double step[3] = { proposed[0] - prevPt[0], proposed[1] - prevPt[1],
          proposed[2] - prevPt[2] };
        double stepLen = vtkMath::Norm(step);
        if (targetEdgeLength > 0.0)
        {
          double maxStep = clampRatio * targetEdgeLength;
          if (stepLen > maxStep && maxStep > 0.0)
          {
            double scale = maxStep / stepLen;
            proposed[0] = prevPt[0] + step[0] * scale;
            proposed[1] = prevPt[1] + step[1] * scale;
            proposed[2] = prevPt[2] + step[2] * scale;
          }
        }
      }

      if (params.source)
      {
        vtkSmoothPoint* sPtr = params.SmoothPoints->GetSmoothPoint(pid);
        vtkCell* cell = nullptr;
        if (sPtr->cellId >= 0)
        {
          cell = params.source->GetCell(sPtr->cellId);
        }

        double closestPt[3];
        double dist2;
        if (!cell ||
          cell->EvaluatePosition(proposed, closestPt, sPtr->subId, sPtr->p, dist2, params.w) == 0)
        {
          params.cellLocator->FindClosestPoint(proposed, closestPt, sPtr->cellId, sPtr->subId, dist2);
        }
        proposed[0] = closestPt[0];
        proposed[1] = closestPt[1];
        proposed[2] = closestPt[2];
      }

      params.newPts->SetPoint(pid, proposed);

      double dispVec[3] = { proposed[0] - prevPt[0], proposed[1] - prevPt[1],
        proposed[2] - prevPt[2] };
      double disp = vtkMath::Norm(dispVec);
      if (disp > maxDist)
      {
        maxDist = disp;
      }
    }

    ++iterationNumber;
    copyCurrentPoints();
  }

  vtkDebugWithObjectMacro(params.spdf,
    << "Performed " << iterationNumber << " uniform triangle smoothing passes");
}

} // namespace

//------------------------------------------------------------------------------
int vtkSmoothPolyDataFilter::RequestData(vtkInformation* vtkNotUsed(request),
  vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
  // get the info objects
  vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation* sourceInfo = inputVector[1]->GetInformationObject(0);
  vtkInformation* outInfo = outputVector->GetInformationObject(0);

  // get the input and output
  vtkPolyData* input = vtkPolyData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkPolyData* source = nullptr;
  if (sourceInfo)
  {
    source = vtkPolyData::SafeDownCast(sourceInfo->Get(vtkDataObject::DATA_OBJECT()));
  }
  vtkPolyData* output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkIdType numPts, numCells, i, numPolys, numStrips;
  int j, k;
  vtkIdType npts = 0;
  const vtkIdType* pts = nullptr;
  vtkIdType p1, p2;
  double conv;
  double x1[3], x2[3], x3[3], l1[3], l2[3];
  double CosFeatureAngle; // Cosine of angle between adjacent polys
  double CosEdgeAngle;    // Cosine of angle between adjacent edges
  double closestPt[3], dist2;
  vtkIdType numSimple = 0, numBEdges = 0, numFixed = 0, numFEdges = 0;
  bool uniformModeActive = (this->SmoothingMode == VTK_SMOOTH_MODE_UNIFORM_TRIANGLE);
  std::vector<std::array<vtkIdType, 3>> uniformTriangles;
  std::vector<std::vector<vtkIdType>> uniformIncidentTriangles;
  vtkPolyData* Mesh;
  vtkPoints* inPts;
  vtkCellArray *inVerts, *inLines, *inPolys, *inStrips;

  // Check input
  //
  numPts = input->GetNumberOfPoints();
  numCells = input->GetNumberOfCells();
  if (numPts < 1 || numCells < 1)
  {
    vtkErrorMacro(<< "No data to smooth!");
    return 1;
  }

  CosFeatureAngle = cos(vtkMath::RadiansFromDegrees(this->FeatureAngle));
  CosEdgeAngle = cos(vtkMath::RadiansFromDegrees(this->EdgeAngle));

  vtkDebugMacro(<< "Smoothing " << numPts << " vertices, " << numCells << " cells with:\n"
                << "\tConvergence= " << this->Convergence << "\n"
                << "\tIterations= " << this->NumberOfIterations << "\n"
                << "\tRelaxation Factor= " << this->RelaxationFactor << "\n"
                << "\tEdge Angle= " << this->EdgeAngle << "\n"
                << "\tBoundary Smoothing " << (this->BoundarySmoothing ? "On\n" : "Off\n")
                << "\tFeature Edge Smoothing " << (this->FeatureEdgeSmoothing ? "On\n" : "Off\n")
                << "\tError Scalars " << (this->GenerateErrorScalars ? "On\n" : "Off\n")
                << "\tError Vectors " << (this->GenerateErrorVectors ? "On\n" : "Off\n"));

  if (this->NumberOfIterations <= 0 || this->RelaxationFactor == 0.0)
  { // don't do anything! pass data through
    output->CopyStructure(input);
    output->GetPointData()->PassData(input->GetPointData());
    output->GetCellData()->PassData(input->GetCellData());
    return 1;
  }

  // Perform topological analysis. What we're gonna do is build a connectivity
  // array of connected vertices. The outcome will be one of three
  // classifications for a vertex: VTK_SIMPLE_VERTEX, VTK_FIXED_VERTEX. or
  // VTK_EDGE_VERTEX. Simple vertices are smoothed using all connected
  // vertices. FIXED vertices are never smoothed. Edge vertices are smoothed
  // using a subset of the attached vertices.
  //
  vtkDebugMacro(<< "Analyzing topology...");

  // Smart pointer to storage; use a raw pointer for operator[] array access.
  std::unique_ptr<vtkMeshVertex[]> uVerts(new vtkMeshVertex[numPts]);
  vtkMeshVertex* Verts = uVerts.get();

  inPts = input->GetPoints();
  conv = this->Convergence * input->GetLength();

  // check vertices first. Vertices are never smoothed_--------------
  for (inVerts = input->GetVerts(), inVerts->InitTraversal(); inVerts->GetNextCell(npts, pts);)
  {
    for (j = 0; j < npts; j++)
    {
      Verts[pts[j]].type = VTK_FIXED_VERTEX;
    }
  }
  this->UpdateProgress(0.10);
  vtkIdType checkAbortInterval = std::min(input->GetNumberOfLines() / 10 + 1, (vtkIdType)1000);
  vtkIdType progressCounter = 0;

  // now check lines. Only manifold lines can be smoothed------------
  for (inLines = input->GetLines(), inLines->InitTraversal(); inLines->GetNextCell(npts, pts);)
  {
    if (progressCounter % checkAbortInterval == 0 && this->CheckAbort())
    {
      break;
    }
    progressCounter++;
    for (j = 0; j < npts; j++)
    {
      if (Verts[pts[j]].type == VTK_SIMPLE_VERTEX)
      {
        if (j == (npts - 1)) // end-of-line marked FIXED
        {
          Verts[pts[j]].type = VTK_FIXED_VERTEX;
        }
        else if (j == 0) // beginning-of-line marked FIXED
        {
          Verts[pts[0]].type = VTK_FIXED_VERTEX;
          inPts->GetPoint(pts[0], x2);
          inPts->GetPoint(pts[1], x3);
        }
        else // is edge vertex (unless already edge vertex!)
        {
          Verts[pts[j]].type = VTK_FEATURE_EDGE_VERTEX;
          Verts[pts[j]].edges = vtkIdList::New();
          Verts[pts[j]].edges->SetNumberOfIds(2);
          Verts[pts[j]].edges->SetId(0, pts[j - 1]);
          Verts[pts[j]].edges->SetId(1, pts[j + 1]);
        }
      } // if simple vertex

      else if (Verts[pts[j]].type == VTK_FEATURE_EDGE_VERTEX)
      { // multiply connected, becomes fixed!
        Verts[pts[j]].type = VTK_FIXED_VERTEX;
        Verts[pts[j]].edges->Delete();
        Verts[pts[j]].edges = nullptr;
      }

    } // for all points in this line
  }   // for all lines
  this->UpdateProgress(0.25);

  // now polygons and triangle strips-------------------------------
  inPolys = input->GetPolys();
  numPolys = inPolys->GetNumberOfCells();
  inStrips = input->GetStrips();
  numStrips = inStrips->GetNumberOfCells();

  if (numPolys > 0 || numStrips > 0)
  { // build cell structure
    vtkCellArray* polys;
    vtkIdType cellId;
    int numNei, nei, edge;
    vtkIdType numNeiPts;
    const vtkIdType* neiPts;
    double normal[3], neiNormal[3];

    vtkNew<vtkIdList> neighbors;
    neighbors->Allocate(VTK_CELL_SIZE);

    vtkNew<vtkPolyData> inMesh;
    inMesh->SetPoints(inPts);
    inMesh->SetPolys(inPolys);
    Mesh = inMesh;

    vtkSmartPointer<vtkTriangleFilter> toTris;
    if ((numStrips = inStrips->GetNumberOfCells()) > 0)
    { // convert data to triangles
      inMesh->SetStrips(inStrips);
      toTris.TakeReference(vtkTriangleFilter::New());
      toTris->SetInputData(inMesh);
      toTris->Update();
      Mesh = toTris->GetOutput();
    }

    Mesh->BuildLinks(); // to do neighborhood searching
    polys = Mesh->GetPolys();
    this->UpdateProgress(0.375);

    checkAbortInterval = std::min(polys->GetNumberOfCells() / 10 + 1, (vtkIdType)1000);

    for (cellId = 0, polys->InitTraversal(); polys->GetNextCell(npts, pts); cellId++)
    {
      if (cellId % checkAbortInterval == 0 && this->CheckAbort())
      {
        break;
      }
      for (i = 0; i < npts; i++)
      {
        p1 = pts[i];
        p2 = pts[(i + 1) % npts];

        if (Verts[p1].edges == nullptr)
        {
          Verts[p1].edges = vtkIdList::New();
          Verts[p1].edges->Allocate(16, 6);
        }
        if (Verts[p2].edges == nullptr)
        {
          Verts[p2].edges = vtkIdList::New();
          Verts[p2].edges->Allocate(16, 6);
        }

        Mesh->GetCellEdgeNeighbors(cellId, p1, p2, neighbors);
        numNei = neighbors->GetNumberOfIds();

        edge = VTK_SIMPLE_VERTEX;
        if (numNei == 0)
        {
          edge = VTK_BOUNDARY_EDGE_VERTEX;
        }

        else if (numNei >= 2)
        {
          // check to make sure that this edge hasn't been marked already
          for (j = 0; j < numNei; j++)
          {
            if (neighbors->GetId(j) < cellId)
            {
              break;
            }
          }
          if (j >= numNei)
          {
            edge = VTK_FEATURE_EDGE_VERTEX;
          }
        }

        else if (numNei == 1 && (nei = neighbors->GetId(0)) > cellId)
        {
          if (this->FeatureEdgeSmoothing)
          {
            vtkPolygon::ComputeNormal(inPts, npts, pts, normal);
            Mesh->GetCellPoints(nei, numNeiPts, neiPts);
            vtkPolygon::ComputeNormal(inPts, numNeiPts, neiPts, neiNormal);

            if (vtkMath::Dot(normal, neiNormal) <= CosFeatureAngle)
            {
              edge = VTK_FEATURE_EDGE_VERTEX;
            }
          }
        }
        else // a visited edge; skip rest of analysis
        {
          continue;
        }

        if (edge && Verts[p1].type == VTK_SIMPLE_VERTEX)
        {
          Verts[p1].edges->Reset();
          Verts[p1].edges->InsertNextId(p2);
          Verts[p1].type = edge;
        }
        else if ((edge && Verts[p1].type == VTK_BOUNDARY_EDGE_VERTEX) ||
          (edge && Verts[p1].type == VTK_FEATURE_EDGE_VERTEX) ||
          (!edge && Verts[p1].type == VTK_SIMPLE_VERTEX))
        {
          Verts[p1].edges->InsertNextId(p2);
          if (Verts[p1].type && edge == VTK_BOUNDARY_EDGE_VERTEX)
          {
            Verts[p1].type = VTK_BOUNDARY_EDGE_VERTEX;
          }
        }

        if (edge && Verts[p2].type == VTK_SIMPLE_VERTEX)
        {
          Verts[p2].edges->Reset();
          Verts[p2].edges->InsertNextId(p1);
          Verts[p2].type = edge;
        }
        else if ((edge && Verts[p2].type == VTK_BOUNDARY_EDGE_VERTEX) ||
          (edge && Verts[p2].type == VTK_FEATURE_EDGE_VERTEX) ||
          (!edge && Verts[p2].type == VTK_SIMPLE_VERTEX))
        {
          Verts[p2].edges->InsertNextId(p1);
          if (Verts[p2].type && edge == VTK_BOUNDARY_EDGE_VERTEX)
          {
            Verts[p2].type = VTK_BOUNDARY_EDGE_VERTEX;
          }
        }
      }
    }
  } // if strips or polys

  this->UpdateProgress(0.50);

  checkAbortInterval = std::min(numPts / 10 + 1, (vtkIdType)1000);

  // post-process edge vertices to make sure we can smooth them
  for (i = 0; i < numPts; i++)
  {
    if (i % checkAbortInterval == 0 && this->CheckAbort())
    {
      break;
    }
    if (Verts[i].type == VTK_SIMPLE_VERTEX)
    {
      numSimple++;
    }

    else if (Verts[i].type == VTK_FIXED_VERTEX)
    {
      numFixed++;
    }

    else if (Verts[i].type == VTK_FEATURE_EDGE_VERTEX || Verts[i].type == VTK_BOUNDARY_EDGE_VERTEX)
    { // see how many edges; if two, what the angle is

      if (!this->BoundarySmoothing && Verts[i].type == VTK_BOUNDARY_EDGE_VERTEX)
      {
        Verts[i].type = VTK_FIXED_VERTEX;
        numBEdges++;
      }

      else if ((npts = Verts[i].edges->GetNumberOfIds()) != 2)
      {
        Verts[i].type = VTK_FIXED_VERTEX;
        numFixed++;
      }

      else // check angle between edges
      {
        inPts->GetPoint(Verts[i].edges->GetId(0), x1);
        inPts->GetPoint(i, x2);
        inPts->GetPoint(Verts[i].edges->GetId(1), x3);

        for (k = 0; k < 3; k++)
        {
          l1[k] = x2[k] - x1[k];
          l2[k] = x3[k] - x2[k];
        }
        if (vtkMath::Normalize(l1) >= 0.0 && vtkMath::Normalize(l2) >= 0.0 &&
          vtkMath::Dot(l1, l2) < CosEdgeAngle)
        {
          numFixed++;
          Verts[i].type = VTK_FIXED_VERTEX;
        }
        else
        {
          if (Verts[i].type == VTK_FEATURE_EDGE_VERTEX)
          {
            numFEdges++;
          }
          else
          {
            numBEdges++;
          }
        }
      } // if along edge
    }   // if edge vertex
  }     // for all points

  vtkDebugMacro(<< "Found\n\t" << numSimple << " simple vertices\n\t" << numFEdges
                << " feature edge vertices\n\t" << numBEdges << " boundary edge vertices\n\t"
                << numFixed << " fixed vertices\n\t");
  (void)numSimple;
  (void)numBEdges;
  (void)numFixed;
  (void)numFEdges;

  if (uniformModeActive)
  {
    uniformIncidentTriangles.clear();
    uniformIncidentTriangles.resize(static_cast<size_t>(numPts));

    auto addTriangle = [&](vtkIdType a, vtkIdType b, vtkIdType c) {
      if (a < 0 || b < 0 || c < 0 || a >= numPts || b >= numPts || c >= numPts)
      {
        return;
      }
      uniformTriangles.push_back({ a, b, c });
      vtkIdType triId = static_cast<vtkIdType>(uniformTriangles.size() - 1);
      uniformIncidentTriangles[static_cast<size_t>(a)].push_back(triId);
      uniformIncidentTriangles[static_cast<size_t>(b)].push_back(triId);
      uniformIncidentTriangles[static_cast<size_t>(c)].push_back(triId);
    };

    vtkCellArray* polysForUniform = input->GetPolys();
    if (polysForUniform)
    {
      polysForUniform->InitTraversal();
      while (polysForUniform->GetNextCell(npts, pts))
      {
        if (npts < 3)
        {
          continue;
        }
        vtkIdType first = pts[0];
        for (vtkIdType idx = 1; idx < npts - 1; ++idx)
        {
          addTriangle(first, pts[idx], pts[idx + 1]);
        }
      }
    }

    vtkCellArray* stripsForUniform = input->GetStrips();
    if (stripsForUniform)
    {
      stripsForUniform->InitTraversal();
      while (stripsForUniform->GetNextCell(npts, pts))
      {
        if (npts < 3)
        {
          continue;
        }
        bool flip = false;
        for (vtkIdType idx = 0; idx < npts - 2; ++idx)
        {
          if (!flip)
          {
            addTriangle(pts[idx], pts[idx + 1], pts[idx + 2]);
          }
          else
          {
            addTriangle(pts[idx + 1], pts[idx], pts[idx + 2]);
          }
          flip = !flip;
        }
      }
    }

    if (uniformTriangles.empty())
    {
      uniformModeActive = false;
      vtkDebugMacro("Uniform triangle smoothing requested but no triangle data was detected;"
                    " falling back to Laplacian mode.");
    }
  }

  vtkDebugMacro(<< "Beginning smoothing iterations...");

  // We've setup the topology...now perform Laplacian smoothing
  //
  vtkNew<vtkPoints> newPts;

  // Set the desired precision for the points in the output.
  if (this->OutputPointsPrecision == vtkAlgorithm::DEFAULT_PRECISION)
  {
    newPts->SetDataType(inPts->GetDataType());
  }
  else if (this->OutputPointsPrecision == vtkAlgorithm::SINGLE_PRECISION)
  {
    newPts->SetDataType(VTK_FLOAT);
  }
  else if (this->OutputPointsPrecision == vtkAlgorithm::DOUBLE_PRECISION)
  {
    newPts->SetDataType(VTK_DOUBLE);
  }

  newPts->SetNumberOfPoints(numPts);

  // If a Source is defined, we do constrained smoothing (that is, points are
  // constrained to the surface of the mesh object).
  std::unique_ptr<double[]> w;
  vtkSmartPointer<vtkCellLocator> cellLocator;
  if (source)
  {
    this->SmoothPoints = std::unique_ptr<vtkSmoothPoints>(new vtkSmoothPoints);
    vtkSmoothPoint* sPtr;
    cellLocator.TakeReference(vtkCellLocator::New());
    auto maxCellSize = source->GetMaxCellSize();
    w.reset(new double[maxCellSize]);
    cellLocator->SetDataSet(source);
    cellLocator->BuildLocator();

    for (i = 0; i < numPts; i++)
    {
      sPtr = this->SmoothPoints->InsertSmoothPoint(i);
      cellLocator->FindClosestPoint(
        inPts->GetPoint(i), closestPt, sPtr->cellId, sPtr->subId, dist2);
      newPts->SetPoint(i, closestPt);
    }
  }
  else // smooth normally
  {
    for (i = 0; i < numPts; i++) // initialize to old coordinates
    {
      newPts->SetPoint(i, inPts->GetPoint(i));
    }
  }

  if (newPts->GetDataType() == VTK_DOUBLE)
  {
    vtkSPDF_InternalParams<double> params = { this, this->NumberOfIterations, newPts,
      this->RelaxationFactor, conv, numPts, Verts, source, this->SmoothPoints.get(), w.get(),
      cellLocator };

    if (uniformModeActive)
    {
      vtkSPDF_MovePointsUniformTriangles(params, uniformTriangles, uniformIncidentTriangles);
    }
    else
    {
      vtkSPDF_MovePoints(params);
    }
  }
  else
  {
    vtkSPDF_InternalParams<float> params = { this, this->NumberOfIterations, newPts,
      static_cast<float>(this->RelaxationFactor), static_cast<float>(conv), numPts, Verts, source,
      this->SmoothPoints.get(), w.get(), cellLocator };

    if (uniformModeActive)
    {
      vtkSPDF_MovePointsUniformTriangles(params, uniformTriangles, uniformIncidentTriangles);
    }
    else
    {
      vtkSPDF_MovePoints(params);
    }
  }

  // Release memory if it's been allocated
  this->SmoothPoints.reset(nullptr);

  // Update output. Only point coordinates have changed.
  //
  output->GetPointData()->PassData(input->GetPointData());
  output->GetCellData()->PassData(input->GetCellData());

  if (this->GenerateErrorScalars)
  {
    vtkNew<vtkFloatArray> newScalars;
    newScalars->SetNumberOfTuples(numPts);
    for (i = 0; i < numPts; i++)
    {
      inPts->GetPoint(i, x1);
      newPts->GetPoint(i, x2);
      newScalars->SetComponent(i, 0, sqrt(vtkMath::Distance2BetweenPoints(x1, x2)));
    }
    int idx = output->GetPointData()->AddArray(newScalars);
    output->GetPointData()->SetActiveAttribute(idx, vtkDataSetAttributes::SCALARS);
  }

  if (this->GenerateErrorVectors)
  {
    vtkNew<vtkFloatArray> newVectors;
    newVectors->SetNumberOfComponents(3);
    newVectors->SetNumberOfTuples(numPts);
    for (i = 0; i < numPts; i++)
    {
      inPts->GetPoint(i, x1);
      newPts->GetPoint(i, x2);
      for (j = 0; j < 3; j++)
      {
        x3[j] = x2[j] - x1[j];
      }
      newVectors->SetTuple(i, x3);
    }
    output->GetPointData()->SetVectors(newVectors);
  }

  output->SetPoints(newPts);

  output->SetVerts(input->GetVerts());
  output->SetLines(input->GetLines());
  output->SetPolys(input->GetPolys());
  output->SetStrips(input->GetStrips());

  // free up connectivity storage
  for (i = 0; i < numPts; i++)
  {
    if (Verts[i].edges)
    {
      Verts[i].edges->Delete();
      Verts[i].edges = nullptr;
    }
  }

  return 1;
}

//------------------------------------------------------------------------------
int vtkSmoothPolyDataFilter::FillInputPortInformation(int port, vtkInformation* info)
{
  if (!this->Superclass::FillInputPortInformation(port, info))
  {
    return 0;
  }

  if (port == 1)
  {
    info->Set(vtkAlgorithm::INPUT_IS_OPTIONAL(), 1);
  }
  return 1;
}

//------------------------------------------------------------------------------
void vtkSmoothPolyDataFilter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);

  os << indent << "Convergence: " << this->Convergence << "\n";
  os << indent << "Number of Iterations: " << this->NumberOfIterations << "\n";
  os << indent << "Relaxation Factor: " << this->RelaxationFactor << "\n";
  os << indent << "Feature Edge Smoothing: " << (this->FeatureEdgeSmoothing ? "On\n" : "Off\n");
  os << indent << "Feature Angle: " << this->FeatureAngle << "\n";
  os << indent << "Edge Angle: " << this->EdgeAngle << "\n";
  os << indent << "Boundary Smoothing: " << (this->BoundarySmoothing ? "On\n" : "Off\n");
  os << indent << "Generate Error Scalars: " << (this->GenerateErrorScalars ? "On\n" : "Off\n");
  os << indent << "Generate Error Vectors: " << (this->GenerateErrorVectors ? "On\n" : "Off\n");
  if (this->GetSource())
  {
    os << indent << "Source: " << static_cast<void*>(this->GetSource()) << "\n";
  }
  else
  {
    os << indent << "Source (none)\n";
  }

  os << indent << "Smoothing Mode: "
     << (this->SmoothingMode == VTK_SMOOTH_MODE_UNIFORM_TRIANGLE ? "UniformTriangle\n"
                                                                 : "Laplacian\n");
  os << indent << "Output Points Precision: " << this->OutputPointsPrecision << "\n";
}
VTK_ABI_NAMESPACE_END
