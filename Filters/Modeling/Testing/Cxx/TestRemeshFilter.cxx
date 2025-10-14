#include "vtkRemeshFilter.h"

#include <vtkNew.h>
#include <vtkOBJReader.h>
#include <vtkPolyData.h>
#include <vtkTestUtilities.h>

#include <cstdlib>

int TestRemeshFilter(int argc, char* argv[])
{
  vtkNew<vtkOBJReader> reader;
  char* fileName = vtkTestUtilities::ExpandDataFileName(argc, argv, "Data/bunny.obj");
  reader->SetFileName(fileName);
  reader->Update();
  delete[] fileName;

  vtkNew<vtkRemeshFilter> remesh;
  remesh->SetInputConnection(reader->GetOutputPort());
  remesh->SetTargetEdgeLength(0.01);
  remesh->SetFeatureAngle(30.0);
  remesh->SetIterations(3);
  remesh->SetProjectToSurface(true);
  remesh->Update();

  vtkPolyData* output = remesh->GetOutput();
  if (!output)
  {
    return EXIT_FAILURE;
  }

  if (output->GetNumberOfPoints() == 0 || output->GetNumberOfCells() == 0)
  {
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
