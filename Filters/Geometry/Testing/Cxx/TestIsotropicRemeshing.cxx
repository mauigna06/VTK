#include "vtkIsotropicRemeshing.h"

#include "vtkCellArray.h"
#include "vtkMath.h"
#include "vtkNew.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkSphereSource.h"
#include "vtkTriangleFilter.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <unordered_map>
#include <vector>

namespace
{
struct EdgeKey
{
  vtkIdType A;
  vtkIdType B;

  bool operator==(const EdgeKey& other) const noexcept
  {
    return this->A == other.A && this->B == other.B;
  }
};

struct EdgeKeyHash
{
  std::size_t operator()(const EdgeKey& k) const noexcept
  {
    std::size_t h1 = std::hash<vtkIdType>{}(k.A);
    std::size_t h2 = std::hash<vtkIdType>{}(k.B);
    return h1 ^ (h2 + 0x9e3779b97f4a7c15ULL + (h1 << 6) + (h1 >> 2));
  }
};

bool ComputeEdgeStats(vtkPolyData* mesh, double& mean, double& stddev)
{
  if (!mesh)
  {
    return false;
  }

  vtkPoints* pts = mesh->GetPoints();
  vtkCellArray* polys = mesh->GetPolys();
  if (!pts || !polys)
  {
    return false;
  }

  std::unordered_map<EdgeKey, double, EdgeKeyHash> edges;
  vtkIdType npts = 0;
  const vtkIdType* ids = nullptr;
  polys->InitTraversal();
  while (polys->GetNextCell(npts, ids))
  {
    if (npts != 3)
    {
      continue;
    }

    for (vtkIdType i = 0; i < 3; ++i)
    {
      vtkIdType a = ids[i];
      vtkIdType b = ids[(i + 1) % 3];
      if (a == b)
      {
        continue;
      }
      EdgeKey key{ std::min(a, b), std::max(a, b) };
      if (edges.find(key) != edges.end())
      {
        continue;
      }

      double pa[3];
      double pb[3];
      pts->GetPoint(a, pa);
      pts->GetPoint(b, pb);
      double length = std::sqrt(vtkMath::Distance2BetweenPoints(pa, pb));
      edges.emplace(key, length);
    }
  }

  if (edges.empty())
  {
    return false;
  }

  std::vector<double> lengths;
  lengths.reserve(edges.size());
  for (const auto& kv : edges)
  {
    lengths.push_back(kv.second);
  }

  mean = std::accumulate(lengths.begin(), lengths.end(), 0.0) /
    static_cast<double>(lengths.size());
  double sumSq = 0.0;
  for (double len : lengths)
  {
    double diff = len - mean;
    sumSq += diff * diff;
  }
  stddev = std::sqrt(sumSq / static_cast<double>(lengths.size()));
  return true;
}

bool ComputeVertexAreaSpread(vtkPolyData* mesh, double& meanArea, double& stdArea)
{
  if (!mesh)
  {
    return false;
  }
  vtkPoints* pts = mesh->GetPoints();
  vtkCellArray* polys = mesh->GetPolys();
  if (!pts || !polys)
  {
    return false;
  }

  std::vector<double> areaPerVertex(static_cast<std::size_t>(pts->GetNumberOfPoints()), 0.0);
  double totalArea = 0.0;
  vtkIdType npts = 0;
  const vtkIdType* ids = nullptr;
  polys->InitTraversal();
  while (polys->GetNextCell(npts, ids))
  {
    if (npts != 3)
    {
      continue;
    }

    double p0[3], p1[3], p2[3];
    pts->GetPoint(ids[0], p0);
    pts->GetPoint(ids[1], p1);
    pts->GetPoint(ids[2], p2);
    double v1[3];
    double v2[3];
    vtkMath::Subtract(p1, p0, v1);
    vtkMath::Subtract(p2, p0, v2);
    double triNormal[3];
    vtkMath::Cross(v1, v2, triNormal);
    double area = 0.5 * vtkMath::Norm(triNormal);
    if (area <= 0.0)
    {
      continue;
    }
    totalArea += area;
    double share = area / 3.0;
    for (vtkIdType i = 0; i < 3; ++i)
    {
      areaPerVertex[ids[i]] += share;
    }
  }

  if (totalArea <= 0.0)
  {
    return false;
  }

  meanArea = totalArea / static_cast<double>(areaPerVertex.size());
  double sumSq = 0.0;
  for (double area : areaPerVertex)
  {
    double diff = area - meanArea;
    sumSq += diff * diff;
  }
  stdArea = std::sqrt(sumSq / static_cast<double>(areaPerVertex.size()));
  return true;
}
}

int TestIsotropicRemeshing(int vtkNotUsed(argc), char* vtkNotUsed(argv)[])
{
  vtkNew<vtkSphereSource> sphere;
  sphere->SetRadius(12.0);
  sphere->SetThetaResolution(24);
  sphere->SetPhiResolution(16);
  sphere->Update();

  vtkNew<vtkTriangleFilter> triangulate;
  triangulate->SetInputConnection(sphere->GetOutputPort());
  triangulate->Update();

  vtkPolyData* input = triangulate->GetOutput();

  double edgeMeanBefore = 0.0;
  double edgeStdBefore = 0.0;
  if (!ComputeEdgeStats(input, edgeMeanBefore, edgeStdBefore))
  {
    std::cerr << "Failed to compute edge statistics on input." << std::endl;
    return EXIT_FAILURE;
  }

  double areaMeanBefore = 0.0;
  double areaStdBefore = 0.0;
  if (!ComputeVertexAreaSpread(input, areaMeanBefore, areaStdBefore))
  {
    std::cerr << "Failed to compute vertex area stats on input." << std::endl;
    return EXIT_FAILURE;
  }

  vtkNew<vtkIsotropicRemeshing> remeshShort;
  remeshShort->SetInputData(input);
  remeshShort->SetNumberOfIterations(25);
  remeshShort->SetRelaxationFactor(0.85);
  remeshShort->SetPreserveTopology(true);
  remeshShort->SetTargetEdgeLength(-1.0);
  remeshShort->Update();

  vtkPolyData* outputShort = remeshShort->GetOutput();

  double edgeMeanShort = 0.0;
  double edgeStdShort = 0.0;
  if (!ComputeEdgeStats(outputShort, edgeMeanShort, edgeStdShort))
  {
    std::cerr << "Failed to compute edge statistics for short run." << std::endl;
    return EXIT_FAILURE;
  }

  double areaMeanShort = 0.0;
  double areaStdShort = 0.0;
  if (!ComputeVertexAreaSpread(outputShort, areaMeanShort, areaStdShort))
  {
    std::cerr << "Failed to compute vertex area stats for short run." << std::endl;
    return EXIT_FAILURE;
  }

  vtkNew<vtkIsotropicRemeshing> remeshLong;
  remeshLong->SetInputData(input);
  remeshLong->SetNumberOfIterations(80);
  remeshLong->SetRelaxationFactor(0.85);
  remeshLong->SetPreserveTopology(true);
  remeshLong->SetTargetEdgeLength(-1.0);
  remeshLong->Update();

  vtkPolyData* outputLong = remeshLong->GetOutput();

  double edgeMeanLong = 0.0;
  double edgeStdLong = 0.0;
  if (!ComputeEdgeStats(outputLong, edgeMeanLong, edgeStdLong))
  {
    std::cerr << "Failed to compute edge statistics on long run output." << std::endl;
    return EXIT_FAILURE;
  }

  double areaMeanLong = 0.0;
  double areaStdLong = 0.0;
  if (!ComputeVertexAreaSpread(outputLong, areaMeanLong, areaStdLong))
  {
    std::cerr << "Failed to compute vertex area stats on long run output." << std::endl;
    return EXIT_FAILURE;
  }

  if (edgeStdShort >= edgeStdBefore)
  {
    std::cerr << "Edge-length variance did not decrease for short run." << std::endl;
    return EXIT_FAILURE;
  }

  if (areaStdShort >= areaStdBefore)
  {
    std::cerr << "Vertex area variance did not decrease for short run." << std::endl;
    return EXIT_FAILURE;
  }

  if (edgeStdLong >= edgeStdShort - 1e-6)
  {
    std::cerr << "Long run did not improve edge variance sufficiently." << std::endl;
    return EXIT_FAILURE;
  }

  if (areaStdLong >= areaStdShort - 1e-6)
  {
    std::cerr << "Long run did not improve vertex area variance sufficiently." << std::endl;
    return EXIT_FAILURE;
  }

  if (std::abs(edgeMeanLong - edgeMeanBefore) / edgeMeanBefore > 0.05)
  {
    std::cerr << "Mean edge length drifted too far from input shape." << std::endl;
    return EXIT_FAILURE;
  }

  if (std::abs(edgeMeanShort - edgeMeanBefore) / edgeMeanBefore > 0.05)
  {
    std::cerr << "Mean edge length drifted too far during short run." << std::endl;
    return EXIT_FAILURE;
  }

  double relEdgeStd = edgeStdLong / edgeMeanLong;
  if (relEdgeStd > 5e-3)
  {
    std::cerr << "Relative edge-length stddev too high: " << relEdgeStd << std::endl;
    return EXIT_FAILURE;
  }

  double relAreaStd = areaStdLong / areaMeanLong;
  if (relAreaStd > 0.08)
  {
    std::cerr << "Relative vertex-area stddev too high: " << relAreaStd << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
