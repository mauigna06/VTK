// SPDX-FileCopyrightText: Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
// SPDX-License-Identifier: BSD-3-Clause
#include "vtkIncenterSubdivisionFilter.h"

#include "vtkCellArray.h"
#include "vtkEdgeTable.h"
#include "vtkIdList.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkCellData.h"
#include "vtkMath.h"
#include <array>
#include <cmath>
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkPointData.h"
#include "vtkSubdivisionFilter.h"

VTK_ABI_NAMESPACE_BEGIN
vtkStandardNewMacro(vtkIncenterSubdivisionFilter);

// We override RequestData to use our custom GenerateSubdivisionCells implementation.
int vtkIncenterSubdivisionFilter::RequestData(
  vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
  // This implementation mirrors vtkInterpolatingSubdivisionFilter::RequestData
  if (!this->Superclass::RequestData(request, inputVector, outputVector))
  {
    return 0;
  }

  // get the info objects
  vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation* outInfo = outputVector->GetInformationObject(0);

  // get the input and output
  vtkPolyData* input = vtkPolyData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkPolyData* output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkIdType numCells;
  int level;
  vtkPoints* outputPts;
  vtkCellArray* outputPolys;
  vtkPointData* outputPD;
  vtkCellData* outputCD;
  vtkIntArray* edgeData;

  vtkPolyData* inputDS = vtkPolyData::New();
  inputDS->CopyStructure(input);
  inputDS->GetPointData()->PassData(input->GetPointData());
  inputDS->GetCellData()->PassData(input->GetCellData());

  for (level = 0; level < this->NumberOfSubdivisions; level++)
  {
    if (this->CheckAbort())
    {
      break;
    }
    inputDS->BuildLinks();
    numCells = inputDS->GetNumberOfCells();

    outputPts = vtkPoints::New();
    outputPts->DeepCopy(inputDS->GetPoints());

    outputPD = vtkPointData::New();
    outputPD->CopyAllocate(inputDS->GetPointData(), 2 * inputDS->GetNumberOfPoints());

    outputCD = vtkCellData::New();
    outputCD->CopyAllocate(inputDS->GetCellData(), 4 * numCells);

    outputPolys = vtkCellArray::New();
    outputPolys->AllocateEstimate(4 * numCells, 3);

    edgeData = vtkIntArray::New();
    edgeData->SetNumberOfComponents(3);
    edgeData->SetNumberOfTuples(numCells);

    if (this->GenerateSubdivisionPoints(inputDS, edgeData, outputPts, outputPD) == 0)
    {
      outputPts->Delete();
      outputPD->Delete();
      outputCD->Delete();
      outputPolys->Delete();
      inputDS->Delete();
      edgeData->Delete();
      vtkErrorMacro("Subdivision failed.");
      return 0;
    }
    // Call our class-specific cell generator
    this->GenerateSubdivisionCells(inputDS, edgeData, outputPolys, outputCD);

    edgeData->Delete();
    inputDS->Delete();
    inputDS = vtkPolyData::New();
    inputDS->SetPoints(outputPts);
    outputPts->Delete();
    inputDS->SetPolys(outputPolys);
    outputPolys->Delete();
    inputDS->GetPointData()->PassData(outputPD);
    outputPD->Delete();
    inputDS->GetCellData()->PassData(outputCD);
    outputCD->Delete();
    inputDS->Squeeze();
  }

  output->SetPoints(inputDS->GetPoints());
  output->SetPolys(inputDS->GetPolys());
  output->GetPointData()->PassData(inputDS->GetPointData());
  output->GetCellData()->PassData(inputDS->GetCellData());
  inputDS->Delete();

  return 1;
}

void vtkIncenterSubdivisionFilter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

int vtkIncenterSubdivisionFilter::GenerateSubdivisionPoints(
  vtkPolyData* inputDS, vtkIntArray* edgeData, vtkPoints* outputPts, vtkPointData* outputPD)
{
  const vtkIdType* pts = nullptr;
  vtkIdType npts, cellId;
  vtkCellArray* inputPolys = inputDS->GetPolys();
  vtkPoints* inputPts = inputDS->GetPoints();
  vtkPointData* inputPD = inputDS->GetPointData();

  double total = inputPolys->GetNumberOfCells();
  double curr = 0;
  bool abort = false;

  // Helper: compute triangle incenter using user-provided logic
  auto compute_triangle_incenter = [](const std::array<double, 3>& A,
                                      const std::array<double, 3>& B,
                                      const std::array<double, 3>& C) {
    auto distance = [&](const std::array<double, 3>& p1, const std::array<double, 3>& p2) {
      return std::sqrt(std::pow(p2[0] - p1[0], 2) + std::pow(p2[1] - p1[1], 2) +
                       std::pow(p2[2] - p1[2], 2));
    };
    double a = distance(B, C);
    double b = distance(C, A);
    double c = distance(A, B);
    double wsum = a + b + c;
    if (wsum == 0.0)
    {
      return std::array<double, 3>{(A[0] + B[0] + C[0]) / 3.0, (A[1] + B[1] + C[1]) / 3.0,
                                    (A[2] + B[2] + C[2]) / 3.0};
    }
    return std::array<double, 3>{(a * A[0] + b * B[0] + c * C[0]) / wsum,
                                 (a * A[1] + b * B[1] + c * C[1]) / wsum,
                                 (a * A[2] + b * B[2] + c * C[2]) / wsum};
  };

  // For each input triangle, compute its incenter, insert it, copy/interpolate
  for (cellId = 0, inputPolys->InitTraversal(); !abort && inputPolys->GetNextCell(npts, pts);
       cellId++)
  {
    if (this->CheckAbort())
    {
      abort = true;
      break;
    }
    if (npts != 3)
    {
      // only triangles supported
      edgeData->SetComponent(cellId, 0, -1);
      edgeData->SetComponent(cellId, 1, -1);
      edgeData->SetComponent(cellId, 2, -1);
      continue;
    }

    double pA[3], pB[3], pC[3];
    inputPts->GetPoint(pts[0], pA);
    inputPts->GetPoint(pts[1], pB);
    inputPts->GetPoint(pts[2], pC);

    std::array<double, 3> A{pA[0], pA[1], pA[2]};
    std::array<double, 3> B{pB[0], pB[1], pB[2]};
    std::array<double, 3> C{pC[0], pC[1], pC[2]};

    auto I = compute_triangle_incenter(A, B, C);

    // Insert the incenter into outputPts
    vtkIdType incenterId = outputPts->InsertNextPoint(I.data());

    // Interpolate point data from the three vertices using incenter barycentric weights
    // weights are proportional to opposite side lengths: a = |BC|, b = |CA|, c = |AB|
    auto distance = [](const std::array<double, 3>& p1, const std::array<double, 3>& p2) {
      return std::sqrt(std::pow(p2[0] - p1[0], 2) + std::pow(p2[1] - p1[1], 2) +
                       std::pow(p2[2] - p1[2], 2));
    };
    double a_len = distance(B, C);
    double b_len = distance(C, A);
    double c_len = distance(A, B);
    double wsum = a_len + b_len + c_len;
    double weights3[3];
    if (wsum == 0.0)
    {
      weights3[0] = weights3[1] = weights3[2] = 1.0 / 3.0;
    }
    else
    {
      weights3[0] = a_len / wsum; // weight for A
      weights3[1] = b_len / wsum; // weight for B
      weights3[2] = c_len / wsum; // weight for C
    }

    vtkSmartPointer<vtkIdList> stencil = vtkSmartPointer<vtkIdList>::New();
    stencil->SetNumberOfIds(3);
    stencil->SetId(0, pts[0]);
    stencil->SetId(1, pts[1]);
    stencil->SetId(2, pts[2]);
    outputPD->InterpolatePoint(inputPD, incenterId, stencil, weights3);

    // Store the incenter index in edgeData so our GenerateSubdivisionCells can use it.
    // We'll store the same incenter id in all three components for convenience.
    edgeData->SetComponent(cellId, 0, static_cast<double>(incenterId));
    edgeData->SetComponent(cellId, 1, static_cast<double>(incenterId));
    edgeData->SetComponent(cellId, 2, static_cast<double>(incenterId));

    this->UpdateProgress(curr / total);
    curr += 1;
  }

  return 1;
}

void vtkIncenterSubdivisionFilter::GenerateSubdivisionCells(
  vtkPolyData* inputDS, vtkIntArray* edgeData, vtkCellArray* outputPolys, vtkCellData* outputCD)
{
  vtkIdType numCells = inputDS->GetNumberOfCells();
  vtkIdType cellId, newId;
  vtkIdType npts;
  const vtkIdType* pts;
  double edgePts[3];
  vtkIdType newCellPts[3];
  vtkCellData* inputCD = inputDS->GetCellData();

  for (cellId = 0; cellId < numCells; cellId++)
  {
    if (inputDS->GetCellType(cellId) != VTK_TRIANGLE)
    {
      continue;
    }
    inputDS->GetCellPoints(cellId, npts, pts);
    edgeData->GetTuple(cellId, edgePts);
    vtkIdType incenterId = static_cast<vtkIdType>(edgePts[0]);

    // Create three triangles: (v0, v1, I), (v1, v2, I), (v2, v0, I)
    newCellPts[0] = pts[0];
    newCellPts[1] = pts[1];
    newCellPts[2] = incenterId;
    newId = outputPolys->InsertNextCell(3, newCellPts);
    outputCD->CopyData(inputCD, cellId, newId);

    newCellPts[0] = pts[1];
    newCellPts[1] = pts[2];
    newCellPts[2] = incenterId;
    newId = outputPolys->InsertNextCell(3, newCellPts);
    outputCD->CopyData(inputCD, cellId, newId);

    newCellPts[0] = pts[2];
    newCellPts[1] = pts[0];
    newCellPts[2] = incenterId;
    newId = outputPolys->InsertNextCell(3, newCellPts);
    outputCD->CopyData(inputCD, cellId, newId);
  }
}
VTK_ABI_NAMESPACE_END
