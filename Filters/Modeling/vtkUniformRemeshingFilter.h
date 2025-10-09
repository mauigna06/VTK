/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkUniformRemeshingFilter.h

=========================================================================*/
#ifndef vtkUniformRemeshingFilter_h
#define vtkUniformRemeshingFilter_h

#include "vtkPolyDataAlgorithm.h"

class vtkIdTypeArray;

class vtkUniformRemeshingFilter : public vtkPolyDataAlgorithm
{
public:
  static vtkUniformRemeshingFilter* New();
  vtkTypeMacro(vtkUniformRemeshingFilter, vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  vtkGetMacro(TargetEdgeLength, double);
  vtkSetMacro(TargetEdgeLength, double);

  vtkGetMacro(MaxIterations, int);
  vtkSetMacro(MaxIterations, int);

protected:
  vtkUniformRemeshingFilter();
  ~vtkUniformRemeshingFilter() override;

  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;

  double TargetEdgeLength;
  int MaxIterations;

private:
  vtkUniformRemeshingFilter(const vtkUniformRemeshingFilter&) = delete;
  void operator=(const vtkUniformRemeshingFilter&) = delete;
};

#endif // vtkUniformRemeshingFilter_h
