// SPDX-FileCopyrightText: Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
// SPDX-License-Identifier: BSD-3-Clause
/**
 * @class   vtkIncenterSubdivisionFilter
 * @brief   generate a subdivision surface using an incenter-based scheme
 *
 * vtkIncenterSubdivisionFilter is a filter that generates output by
 * subdividing its input polydata. Each subdivision iteration creates new
 * points on edges using the incenter of the edge endpoints (here implemented
 * as the midpoint for a simple example — behavior mirrors the linear scheme
 * but class is provided for specialization).
 *
 * @sa
 * vtkInterpolatingSubdivisionFilter vtkLinearSubdivisionFilter
 */

#ifndef vtkIncenterSubdivisionFilter_h
#define vtkIncenterSubdivisionFilter_h

#include "vtkFiltersModelingModule.h" // For export macro
#include "vtkInterpolatingSubdivisionFilter.h"

VTK_ABI_NAMESPACE_BEGIN
class vtkIntArray;
class vtkPointData;
class vtkPoints;
class vtkPolyData;

class VTKFILTERSMODELING_EXPORT vtkIncenterSubdivisionFilter
  : public vtkInterpolatingSubdivisionFilter
{
public:
  ///@{
  /**
   * Construct object with NumberOfSubdivisions set to 1.
   */
  static vtkIncenterSubdivisionFilter* New();
  vtkTypeMacro(vtkIncenterSubdivisionFilter, vtkInterpolatingSubdivisionFilter);
  void PrintSelf(ostream& os, vtkIndent indent) override;
  ///@}

protected:
  vtkIncenterSubdivisionFilter() = default;
  ~vtkIncenterSubdivisionFilter() override = default;

  int GenerateSubdivisionPoints(vtkPolyData* inputDS, vtkIntArray* edgeData, vtkPoints* outputPts,
    vtkPointData* outputPD) override;
  // Note: vtkInterpolatingSubdivisionFilter::GenerateSubdivisionCells is not virtual,
  // so we cannot use 'override' here. We implement a class-specific method and
  // override RequestData to call it.
  void GenerateSubdivisionCells(vtkPolyData* inputDS, vtkIntArray* edgeData, vtkCellArray* outputPolys,
    vtkCellData* outputCD);

  int RequestData(vtkInformation* request, vtkInformationVector** inputVector,
    vtkInformationVector* outputVector) override;

private:
  vtkIncenterSubdivisionFilter(const vtkIncenterSubdivisionFilter&) = delete;
  void operator=(const vtkIncenterSubdivisionFilter&) = delete;
};

VTK_ABI_NAMESPACE_END
#endif
