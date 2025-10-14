set(classes
  vtkAdaptiveSubdivisionFilter
  vtkBandedPolyDataContourFilter
  vtkButterflySubdivisionFilter
  vtkCollisionDetectionFilter
  vtkContourLoopExtraction
  vtkCookieCutter
  vtkDijkstraGraphGeodesicPath
  vtkDijkstraImageGeodesicPath
  vtkFillHolesFilter
  vtkFitToHeightMapFilter
  vtkGeodesicPath
  vtkGraphGeodesicPath
  vtkHausdorffDistancePointSetFilter
  vtkHyperTreeGridOutlineFilter
  vtkImageDataOutlineFilter
  vtkImprintFilter
  vtkLinearCellExtrusionFilter
  vtkLinearExtrusionFilter
  vtkLinearSubdivisionFilter
  vtkLoopSubdivisionFilter
  vtkOutlineFilter
  vtkPolyDataPointSampler
  vtkProjectedTexture
  vtkQuadRotationalExtrusionFilter
  vtkRibbonFilter
  vtkRotationalExtrusionFilter
  vtkRuledSurfaceFilter
  vtkSectorSource
  vtkSelectEnclosedPoints
  vtkSelectPolyData
  vtkSpherePuzzle
  vtkSpherePuzzleArrows
  vtkSubdivideTetra
  vtkTrimmedExtrusionFilter
  vtkVolumeOfRevolutionFilter)

set(remesh_sources
  Remesh/BoundingBox.cpp
  Remesh/Bvh.cpp
  Remesh/Edge.cpp
  Remesh/Face.cpp
  Remesh/HalfEdge.cpp
  Remesh/Mesh.cpp
  Remesh/MeshIO.cpp
  Remesh/Vertex.cpp)

list(APPEND classes vtkRemeshFilter)

vtk_module_add_module(VTK::FiltersModeling
  CLASSES ${classes}
  SOURCES ${remesh_sources})
vtk_add_test_mangling(VTK::FiltersModeling)
