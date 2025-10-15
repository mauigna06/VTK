#include "vtkRemeshFilter.h"

#include "Remesh/Mesh.h"
#include "Remesh/MeshIO.h"

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkFieldData.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkObjectFactory.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkTriangleFilter.h>

#include <vtkNew.h>

#include <algorithm>
#include <cstdint>
#include <cmath>
#include <unordered_set>

VTK_ABI_NAMESPACE_BEGIN
vtkStandardNewMacro(vtkRemeshFilter);

vtkRemeshFilter::vtkRemeshFilter()
  : TargetEdgeLength(0.0)
  , FeatureAngle(30.0)
  , Iterations(5)
  , ProjectToSurface(true)
{
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

void vtkRemeshFilter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
  os << indent << "TargetEdgeLength: " << this->TargetEdgeLength << "\n";
  os << indent << "FeatureAngle: " << this->FeatureAngle << "\n";
  os << indent << "Iterations: " << this->Iterations << "\n";
  os << indent << "ProjectToSurface: " << (this->ProjectToSurface ? "On" : "Off") << "\n";
}

int vtkRemeshFilter::FillInputPortInformation(int port, vtkInformation* info)
{
  if (port == 0)
  {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
    return 1;
  }
  return 0;
}

int vtkRemeshFilter::RequestData(
  vtkInformation*, vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
  vtkPolyData* input = vtkPolyData::GetData(inputVector[0], 0);
  vtkPolyData* output = vtkPolyData::GetData(outputVector, 0);

  if (!input || !output)
  {
    return 0;
  }

  vtkNew<vtkTriangleFilter> triangulator;
  triangulator->SetInputData(input);
  triangulator->PassVertsOff();
  triangulator->PassLinesOff();
  triangulator->Update();

  vtkPolyData* workMesh = triangulator->GetOutput();

  vtkBotschKobbeltRemeshing::Mesh halfEdgeMesh;
  if (!vtkBotschKobbeltRemeshing::MeshIO::BuildFromPolyData(workMesh, halfEdgeMesh))
  {
    vtkWarningMacro(<< "Failed to convert input to remeshing data structure. Returning input.");
    output->ShallowCopy(input);
    return 1;
  }

  double edgeLength = this->TargetEdgeLength;
  if (edgeLength <= 0.0)
  {
    edgeLength = this->ComputeDefaultEdgeLength(workMesh);
  }

  if (edgeLength <= 0.0)
  {
    vtkWarningMacro(<< "Unable to determine target edge length. Returning input.");
    output->ShallowCopy(input);
    return 1;
  }

  halfEdgeMesh.remesh(edgeLength, this->FeatureAngle, this->Iterations, this->ProjectToSurface);

  output->Initialize();
  vtkBotschKobbeltRemeshing::MeshIO::CopyToPolyData(halfEdgeMesh, output);
  output->GetPointData()->Initialize();
  output->GetCellData()->Initialize();
  if (input->GetFieldData())
  {
    output->GetFieldData()->ShallowCopy(input->GetFieldData());
  }

  return 1;
}

double vtkRemeshFilter::ComputeDefaultEdgeLength(vtkPolyData* mesh) const
{
  if (!mesh || !mesh->GetPoints())
  {
    return 0.0;
  }

  vtkPoints* points = mesh->GetPoints();
  vtkCellArray* polys = mesh->GetPolys();
  if (!polys)
  {
    return 0.0;
  }

  auto encodeEdge = [](vtkIdType a, vtkIdType b) -> std::uint64_t {
    if (a > b)
    {
      std::swap(a, b);
    }
    return (static_cast<std::uint64_t>(a) << 32) | static_cast<std::uint64_t>(b);
  };

  std::unordered_set<std::uint64_t> uniqueEdges;
  uniqueEdges.reserve(static_cast<std::size_t>(polys->GetNumberOfCells()) * 3);

  polys->InitTraversal();
  vtkIdType npts = 0;
  const vtkIdType* pts = nullptr;
  double total = 0.0;
  vtkIdType edgeCount = 0;

  double p0[3];
  double p1[3];
  while (polys->GetNextCell(npts, pts))
  {
    if (npts < 2)
    {
      continue;
    }

    for (vtkIdType i = 0; i < npts; ++i)
    {
      vtkIdType id0 = pts[i];
      vtkIdType id1 = pts[(i + 1) % npts];
      const std::uint64_t key = encodeEdge(id0, id1);
      if (uniqueEdges.insert(key).second)
      {
        points->GetPoint(id0, p0);
        points->GetPoint(id1, p1);
        double dx = p0[0] - p1[0];
        double dy = p0[1] - p1[1];
        double dz = p0[2] - p1[2];
        total += std::sqrt(dx * dx + dy * dy + dz * dz);
        edgeCount++;
      }
    }
  }

  return edgeCount > 0 ? total / static_cast<double>(edgeCount) : 0.0;
}

VTK_ABI_NAMESPACE_END
