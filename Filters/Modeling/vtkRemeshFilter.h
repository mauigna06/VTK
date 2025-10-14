// SPDX-FileCopyrightText: Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
// SPDX-License-Identifier: BSD-3-Clause
/**
 * @class   vtkRemeshFilter
 * @brief   surface remeshing filter based on Botsch and Kobbelt
 *
 * vtkRemeshFilter converts the input polygonal surface into the internal
 * half-edge data structures provided by the Botsch and Kobbelt remeshing
 * implementation and executes the adaptive remeshing procedure. The filter can
 * either use a user-provided target edge length or estimate one from the input
 * mesh when the target length is set to zero. Only triangle meshes are
 * supported; non-triangular inputs are internally triangulated.
 */

#ifndef vtkRemeshFilter_h
#define vtkRemeshFilter_h

#include "vtkFiltersModelingModule.h" // For export macro
#include "vtkPolyDataAlgorithm.h"

VTK_ABI_NAMESPACE_BEGIN

class VTKFILTERSMODELING_EXPORT vtkRemeshFilter : public vtkPolyDataAlgorithm
{
public:
  static vtkRemeshFilter* New();
  vtkTypeMacro(vtkRemeshFilter, vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  ///@{
  /**
   * Target edge length in world coordinates. A value <= 0 triggers automatic
   * estimation from the input mesh.
   */
  vtkSetMacro(TargetEdgeLength, double);
  vtkGetMacro(TargetEdgeLength, double);
  ///@}

  ///@{
  /**
   * Feature angle in degrees used to preserve sharp features (0-180).
   */
  vtkSetClampMacro(FeatureAngle, double, 0.0, 180.0);
  vtkGetMacro(FeatureAngle, double);
  ///@}

  ///@{
  /**
   * Number of remeshing iterations to execute (minimum 1).
   */
  vtkSetClampMacro(Iterations, int, 1, VTK_INT_MAX);
  vtkGetMacro(Iterations, int);
  ///@}

  ///@{
  /**
   * Enable projection of the remeshed surface back onto the original input.
   */
  vtkSetMacro(ProjectToSurface, bool);
  vtkGetMacro(ProjectToSurface, bool);
  vtkBooleanMacro(ProjectToSurface, bool);
  ///@}

protected:
  vtkRemeshFilter();
  ~vtkRemeshFilter() override = default;

  int FillInputPortInformation(int port, vtkInformation* info) override;
  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;

  double ComputeDefaultEdgeLength(vtkPolyData* mesh) const;

  double TargetEdgeLength;
  double FeatureAngle;
  int Iterations;
  bool ProjectToSurface;

private:
  vtkRemeshFilter(const vtkRemeshFilter&) = delete;
  void operator=(const vtkRemeshFilter&) = delete;
};

VTK_ABI_NAMESPACE_END

#endif
