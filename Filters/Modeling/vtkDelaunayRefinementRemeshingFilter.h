// SPDX-FileCopyrightText: Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
// SPDX-License-Identifier: BSD-3-Clause
/**
 * @class   vtkDelaunayRefinementRemeshingFilter
 * @brief   remesh a surface using Delaunay triangulation + local refinement
 *
 * vtkDelaunayRefinementRemeshingFilter is a scaffolded filter that performs
 * remeshing by creating a Delaunay triangulation of input points and then
 * applying local refinement (edge splitting/removal) to achieve a target edge
 * length and mesh quality. The implementation here uses VTK's
 * vtkDelaunay2D for triangulation; refinement steps are placeholders you can
 * extend (edge length checks, edge flips, Laplacian smoothing, etc.).
 */

#ifndef vtkDelaunayRefinementRemeshingFilter_h
#define vtkDelaunayRefinementRemeshingFilter_h

#include "vtkFiltersModelingModule.h" // For export macro
#include "vtkPolyDataAlgorithm.h"

VTK_ABI_NAMESPACE_BEGIN

class vtkPolyData;

class VTKFILTERSMODELING_EXPORT vtkDelaunayRefinementRemeshingFilter : public vtkPolyDataAlgorithm
{
public:
  static vtkDelaunayRefinementRemeshingFilter* New();
  vtkTypeMacro(vtkDelaunayRefinementRemeshingFilter, vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  // Target edge length for refinement (approximate)
  vtkSetMacro(TargetEdgeLength, double);
  vtkGetMacro(TargetEdgeLength, double);

  // Minimum allowed triangle angle (degrees). Triangles with smaller angles
  // will be refined (Ruppert-style refinement).
  vtkSetMacro(MinAngleDegrees, double);
  vtkGetMacro(MinAngleDegrees, double);

  // Maximum number of refinement iterations
  vtkSetMacro(MaxIterations, int);
  vtkGetMacro(MaxIterations, int);

protected:
  vtkDelaunayRefinementRemeshingFilter();
  ~vtkDelaunayRefinementRemeshingFilter() override;

  int RequestData(vtkInformation* request, vtkInformationVector** inputVector,
    vtkInformationVector* outputVector) override;

private:
  double TargetEdgeLength;
  int MaxIterations;
  double MinAngleDegrees;

  vtkDelaunayRefinementRemeshingFilter(const vtkDelaunayRefinementRemeshingFilter&) = delete;
  void operator=(const vtkDelaunayRefinementRemeshingFilter&) = delete;
};

VTK_ABI_NAMESPACE_END

#endif
