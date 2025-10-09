vtkIsotropicRemeshing
======================

Overview
--------
`vtkIsotropicRemeshing` is a simple vtkPolyDataAlgorithm that performs
iterative centroidal-Voronoi-like relaxation on triangle meshes. It
approximates Lloyd iterations by moving vertices toward the area-weighted
centroid of incident triangle centroids. This improves triangle shape
isotropy but does not implement topological operations (edge splits/collapses)
or an exact surface Voronoi diagram.

Usage
-----
- Include the header: `#include "vtkIsotropicRemeshing.h"`
- Create and configure the filter like any VTK filter:

  vtkNew<vtkIsotropicRemeshing> remesh;
  remesh->SetNumberOfIterations(10);
  remesh->SetRelaxationFactor(0.8);
  remesh->SetPreserveTopology(true);
  remesh->SetInputData(inputPolyData);
  remesh->Update();

Limitations
-----------
- Does not change mesh connectivity. For full isotropic remeshing (edge
  flipping, split/collapse for target edge length) additional topology
  operations are required.
- Boundary vertices are optionally preserved; this filter does not
  reconstruct sharp creases or features.

Notes
-----
This implementation is intended as a compact example of using centroidal
relaxation inside a VTK filter. It is suitable for small-to-medium meshes
and as a starting point for more advanced remeshing implementations.
