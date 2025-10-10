/*=========================================================================

  Module:    vtkIsotropicRemeshing.h

  A lightweight isotropic remeshing filter for vtkPolyData that performs
  centroidal-Voronoi-like relaxation (Lloyd iterations approximated by
  area-weighted centroids on the one-ring) to improve triangle isotropy and
  gently drift samples from locally dense regions toward sparse regions while
  staying on the original surface.

  NOTE: This implementation is a self-contained approximation of centroidal
  Voronoi tessellation suitable as an example or a starting point. It does
  not implement full edge split/collapse or an exact surface Voronoi.

=========================================================================*/
#ifndef vtkIsotropicRemeshing_h
#define vtkIsotropicRemeshing_h

#include "vtkPolyDataAlgorithm.h"
#include "vtkFiltersGeometryModule.h" // For export macro

class VTKFILTERSGEOMETRY_EXPORT vtkIsotropicRemeshing : public vtkPolyDataAlgorithm
{
public:
  static vtkIsotropicRemeshing* New();
  vtkTypeMacro(vtkIsotropicRemeshing, vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  // Number of smoothing iterations (Lloyd-style)
  vtkSetMacro(NumberOfIterations, int);
  vtkGetMacro(NumberOfIterations, int);

  // Relaxation factor (0..1). 1.0 moves the vertex fully to the centroid.
  vtkSetMacro(RelaxationFactor, double);
  vtkGetMacro(RelaxationFactor, double);

  // If <= 0 the target edge length is estimated from the input (average)
  vtkSetMacro(TargetEdgeLength, double);
  vtkGetMacro(TargetEdgeLength, double);

  // Preserve boundary points (do not move boundary vertices)
  vtkSetMacro(PreserveTopology, bool);
  vtkGetMacro(PreserveTopology, bool);
  vtkBooleanMacro(PreserveTopology, bool);

protected:
  vtkIsotropicRemeshing();
  ~vtkIsotropicRemeshing() override;

  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;

  int NumberOfIterations;
  double RelaxationFactor;
  double TargetEdgeLength;
  bool PreserveTopology;

private:
  vtkIsotropicRemeshing(const vtkIsotropicRemeshing&) = delete;
  void operator=(const vtkIsotropicRemeshing&) = delete;
};

#endif // vtkIsotropicRemeshing_h
