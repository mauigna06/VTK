vtkIsotropicRemeshing
======================

Overview
--------
`vtkIsotropicRemeshing` is a vtkPolyDataAlgorithm that performs
iterative centroidal-Voronoi-like relaxation on triangle meshes. It
approximates Lloyd iterations by moving vertices tangentially toward the
area-weighted centroid of their incident triangles while simultaneously
balancing one-ring edge lengths against a target value. The result is a
noticeable improvement in triangle isotropy without altering mesh
connectivity.

Usage
-----
- Include the header: `#include "vtkIsotropicRemeshing.h"`
- Create and configure the filter like any VTK filter:

  vtkNew<vtkIsotropicRemeshing> remesh;
  remesh->SetNumberOfIterations(10);
  remesh->SetRelaxationFactor(0.8);
  remesh->SetPreserveTopology(true);
  remesh->SetTargetEdgeLength(0.0); // optional: auto-estimate from input
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
This implementation is intended as a compact example of centroidal
relaxation paired with edge-length equalisation inside a VTK filter. It
is suitable for small-to-medium meshes and as a foundation for more
advanced remeshing implementations.
