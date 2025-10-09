// SPDX-FileCopyrightText: Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
// SPDX-License-Identifier: BSD-3-Clause
#include "vtkDelaunayRefinementRemeshingFilter.h"

#include "vtkCellArray.h"
#include "vtkDelaunay2D.h"
#include "vtkIdList.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkPointLocator.h"
#include "vtkMath.h"
#include "vtkCell.h"

#include <cmath>
#include <limits>
#include <vector>

VTK_ABI_NAMESPACE_BEGIN

vtkStandardNewMacro(vtkDelaunayRefinementRemeshingFilter);

vtkDelaunayRefinementRemeshingFilter::vtkDelaunayRefinementRemeshingFilter()
  : TargetEdgeLength(1.0)
  , MaxIterations(10)
  , MinAngleDegrees(20.0)
{
}

vtkDelaunayRefinementRemeshingFilter::~vtkDelaunayRefinementRemeshingFilter() = default;

void vtkDelaunayRefinementRemeshingFilter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
  os << indent << "TargetEdgeLength: " << this->TargetEdgeLength << "\n";
  os << indent << "MaxIterations: " << this->MaxIterations << "\n";
  os << indent << "MinAngleDegrees: " << this->MinAngleDegrees << "\n";
}

int vtkDelaunayRefinementRemeshingFilter::RequestData(vtkInformation* vtkNotUsed(request),
  vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
  vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation* outInfo = outputVector->GetInformationObject(0);

  vtkPolyData* input = vtkPolyData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkPolyData* output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

  if (!input || !output)
  {
    return 0;
  }

  // Copy input points to a working set (we'll iteratively add points)
  vtkSmartPointer<vtkPoints> workingPoints = vtkSmartPointer<vtkPoints>::New();
  vtkIdType numInPts = 0;
  if (input->GetPoints())
  {
    numInPts = input->GetPoints()->GetNumberOfPoints();
    for (vtkIdType i = 0; i < numInPts; ++i)
    {
      double p[3];
      input->GetPoints()->GetPoint(i, p);
      workingPoints->InsertNextPoint(p);
    }
  }

  auto makePolyFromPoints = [&](vtkPoints* pts) -> vtkSmartPointer<vtkPolyData> {
    vtkSmartPointer<vtkPolyData> pd = vtkSmartPointer<vtkPolyData>::New();
    pd->SetPoints(pts);
    return pd;
  };

  const double eps = std::numeric_limits<double>::epsilon() * 100.0;

  bool converged = false;
  for (int iter = 0; iter < this->MaxIterations && !converged; ++iter)
  {
    vtkSmartPointer<vtkPolyData> ptsPoly = makePolyFromPoints(workingPoints);

    // Build a point locator to avoid inserting duplicate/nearby points
    vtkSmartPointer<vtkPointLocator> locator = vtkSmartPointer<vtkPointLocator>::New();
    locator->SetDataSet(ptsPoly);
    locator->BuildLocator();

    // Triangulate current point set
    vtkSmartPointer<vtkDelaunay2D> delaunay = vtkSmartPointer<vtkDelaunay2D>::New();
    delaunay->SetInputData(ptsPoly);
    delaunay->Update();

    vtkSmartPointer<vtkPolyData> tri = delaunay->GetOutput();
    vtkIdType ncells = tri->GetNumberOfCells();
    if (ncells == 0)
    {
      // nothing to refine
      break;
    }

    // Track the worst quality triangle (smallest angle) and the longest edge
    double worstAngle = std::numeric_limits<double>::infinity();
    vtkIdType worstCellId = -1;
    double worstEdgeLen = 0.0;
    vtkIdType worstEdgeP0 = -1, worstEdgeP1 = -1;

    vtkSmartPointer<vtkIdList> pts = vtkSmartPointer<vtkIdList>::New();
    for (vtkIdType cid = 0; cid < ncells; ++cid)
    {
      tri->GetCellPoints(cid, pts);
      if (pts->GetNumberOfIds() != 3)
      {
        continue;
      }
      double p0[3], p1[3], p2[3];
      tri->GetPoint(pts->GetId(0), p0);
      tri->GetPoint(pts->GetId(1), p1);
      tri->GetPoint(pts->GetId(2), p2);

      // edge lengths
      double e0 = vtkMath::Distance2BetweenPoints(p0, p1);
      double e1 = vtkMath::Distance2BetweenPoints(p1, p2);
      double e2 = vtkMath::Distance2BetweenPoints(p2, p0);

      double edges[3] = { e0, e1, e2 };
      double ptsArr[3][3];
      for (int k = 0; k < 3; ++k)
      {
        tri->GetPoint(pts->GetId(k), ptsArr[k]);
      }

      // compute angles at each vertex
      double angles[3];
      for (int k = 0; k < 3; ++k)
      {
        double va[3], vb[3];
        vtkMath::Subtract(ptsArr[(k + 1) % 3], ptsArr[k], va);
        vtkMath::Subtract(ptsArr[(k + 2) % 3], ptsArr[k], vb);
        vtkMath::Normalize(va);
        vtkMath::Normalize(vb);
        double dot = vtkMath::Dot(va, vb);
        dot = std::min(1.0, std::max(-1.0, dot));
        angles[k] = vtkMath::DegreesFromRadians(acos(dot));
      }

      double localMinAngle = std::min(std::min(angles[0], angles[1]), angles[2]);
      if (localMinAngle < worstAngle)
      {
        worstAngle = localMinAngle;
        worstCellId = cid;
      }

      // check longest edge (in world length)
      double maxEdge = sqrt(std::max(std::max(e0, e1), e2));
      if (maxEdge > worstEdgeLen)
      {
        worstEdgeLen = maxEdge;
        // store endpoints of the longest edge
        if (e0 >= e1 && e0 >= e2)
        {
          worstEdgeP0 = pts->GetId(0);
          worstEdgeP1 = pts->GetId(1);
        }
        else if (e1 >= e0 && e1 >= e2)
        {
          worstEdgeP0 = pts->GetId(1);
          worstEdgeP1 = pts->GetId(2);
        }
        else
        {
          worstEdgeP0 = pts->GetId(2);
          worstEdgeP1 = pts->GetId(0);
        }
      }
    }

    // Decide on insertion: first fix small angles (Ruppert), else split long edges
    double minAngleThreshold = this->MinAngleDegrees;
    bool inserted = false;
    double insertPt[3];

    if (worstCellId >= 0 && worstAngle < minAngleThreshold)
    {
      // compute circumcenter of the worst triangle
      tri->GetCellPoints(worstCellId, pts);
      double a[3], b[3], cpt[3];
      tri->GetPoint(pts->GetId(0), a);
      tri->GetPoint(pts->GetId(1), b);
      tri->GetPoint(pts->GetId(2), cpt);

      // build local 2D frame on triangle plane
      double ab[3], ac[3];
      vtkMath::Subtract(b, a, ab);
      vtkMath::Subtract(cpt, a, ac);
      double d = vtkMath::Norm(ab);
      if (d < eps)
      {
        // degenerate, skip
      }
      else
      {
        double u[3];
        vtkMath::Normalize(ab);
        u[0] = ab[0];
        u[1] = ab[1];
        u[2] = ab[2];
        double proj = vtkMath::Dot(ac, u);
        double vtemp[3];
        for (int ii = 0; ii < 3; ++ii)
        {
          vtemp[ii] = ac[ii] - proj * u[ii];
        }
        double vnorm = vtkMath::Norm(vtemp);
        if (vnorm < eps)
        {
          // nearly colinear, skip
        }
        else
        {
          double v[3];
          for (int ii = 0; ii < 3; ++ii)
          {
            v[ii] = vtemp[ii] / vnorm;
          }

          // 2D coordinates
          double ax = 0.0, ay = 0.0;
          double bx = d, by = 0.0;
          double cx = proj, cy = vnorm;

          double x1 = ax, y1 = ay;
          double x2 = bx, y2 = by;
          double x3 = cx, y3 = cy;
          double denom = 2.0 * (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2));
          if (std::abs(denom) < eps)
          {
            // degenerate
          }
          else
          {
            double x1s = x1 * x1 + y1 * y1;
            double x2s = x2 * x2 + y2 * y2;
            double x3s = x3 * x3 + y3 * y3;
            double ux = (x1s * (y2 - y3) + x2s * (y3 - y1) + x3s * (y1 - y2)) / denom;
            double uy = (x1s * (x3 - x2) + x2s * (x1 - x3) + x3s * (x2 - x1)) / denom;

            // map back to 3D: a + ux * u + uy * v
            for (int ii = 0; ii < 3; ++ii)
            {
              insertPt[ii] = a[ii] + ux * u[ii] + uy * v[ii];
            }
            inserted = true;
          }
        }
      }
    }
    else if (worstEdgeLen > this->TargetEdgeLength)
    {
      // split the longest edge
      vtkIdType pid0 = worstEdgeP0;
      vtkIdType pid1 = worstEdgeP1;
      if (pid0 >= 0 && pid1 >= 0)
      {
        double p0[3], p1[3];
        tri->GetPoint(pid0, p0);
        tri->GetPoint(pid1, p1);
        for (int ii = 0; ii < 3; ++ii)
        {
          insertPt[ii] = 0.5 * (p0[ii] + p1[ii]);
        }
        inserted = true;
      }
    }

    if (!inserted)
    {
      converged = true;
      break;
    }

    // check distance to nearest existing point
    double closestPt[3];
    vtkIdType nearestId = locator->FindClosestPoint(insertPt);
    // vtkPointLocator does not provide GetClosestPoint; use point id to fetch coords
    if (nearestId >= 0)
    {
      ptsPoly->GetPoint(nearestId, closestPt);
    }
    else
    {
      closestPt[0] = closestPt[1] = closestPt[2] = 0.0;
    }
    double d2 = vtkMath::Distance2BetweenPoints(closestPt, insertPt);
    double minDist = (this->TargetEdgeLength > 0.0) ? (this->TargetEdgeLength * 1e-6) : 1e-8;
    if (d2 <= (minDist * minDist))
    {
      // too close to existing point; abort this insertion to avoid infinite loop
      converged = true;
      break;
    }

    // actually insert
    workingPoints->InsertNextPoint(insertPt);
  }

  // Final triangulation
  vtkSmartPointer<vtkPolyData> finalPts = vtkSmartPointer<vtkPolyData>::New();
  finalPts->SetPoints(workingPoints);
  vtkSmartPointer<vtkDelaunay2D> finalD = vtkSmartPointer<vtkDelaunay2D>::New();
  finalD->SetInputData(finalPts);
  finalD->Update();

  vtkSmartPointer<vtkPolyData> finalTri = finalD->GetOutput();
  output->ShallowCopy(finalTri);

  return 1;
}

VTK_ABI_NAMESPACE_END
