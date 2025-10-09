#include "vtkUniformRemeshingFilter.h"

#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkPointData.h"
#include "vtkCellArray.h"
#include "vtkTriangleFilter.h"
#include "vtkCellLocator.h"
#include "vtkIdList.h"
#include "vtkMath.h"

#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <cmath>
#include <algorithm>

vtkStandardNewMacro(vtkUniformRemeshingFilter);

vtkUniformRemeshingFilter::vtkUniformRemeshingFilter()
{
  this->TargetEdgeLength = 0.0; // disabled by default
  this->MaxIterations = 10;
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

vtkUniformRemeshingFilter::~vtkUniformRemeshingFilter() = default;

void vtkUniformRemeshingFilter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
  os << indent << "TargetEdgeLength: " << this->TargetEdgeLength << "\n";
  os << indent << "MaxIterations: " << this->MaxIterations << "\n";
}

int vtkUniformRemeshingFilter::RequestData(vtkInformation* vtkNotUsed(request),
                                            vtkInformationVector** inputVector,
                                            vtkInformationVector* outputVector)
{
  vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation* outInfo = outputVector->GetInformationObject(0);

  vtkPolyData* input = vtkPolyData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkPolyData* output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

  if (!input)
  {
    vtkErrorMacro("No input polydata");
    return 0;
  }

  // If TargetEdgeLength <= 0, compute average edge length and use it as target
  double target = this->TargetEdgeLength;

  // Triangulate the input to have consistent triangles for remeshing
  vtkNew<vtkTriangleFilter> triangulator;
  triangulator->SetInputData(input);
  triangulator->PassVertsOff();
  triangulator->PassLinesOff();
  triangulator->Update();

  vtkSmartPointer<vtkPolyData> working = vtkSmartPointer<vtkPolyData>::New();
  working->DeepCopy(triangulator->GetOutput());

  if (working->GetNumberOfCells() == 0 || working->GetNumberOfPoints() == 0)
  {
    vtkErrorMacro("Input contains no geometry after triangulation");
    return 0;
  }

  // Non-CGAL VTK-only isotropic-style remeshing
  // Build adjacency: for each point, store neighboring point ids and incident cells
  vtkIdType np = working->GetNumberOfPoints();
  vtkIdType nc = working->GetNumberOfCells();

  // Build edge list and adjacency
  struct Edge { vtkIdType a, b; double len; };
  std::vector<Edge> edges;
  edges.reserve(working->GetNumberOfCells()*3);

  vtkCellArray* polys = working->GetPolys();
  polys->InitTraversal();
  vtkIdType n; const vtkIdType* ids;
  std::unordered_map<uint64_t, vtkIdType> edgeMap; // map ordered pair to edge index
  auto edgeKey = [](vtkIdType x, vtkIdType y)->uint64_t { return (static_cast<uint64_t>(std::min(x,y))<<32) | static_cast<uint64_t>(std::max(x,y)); };

  vtkPoints* pts = working->GetPoints();
  while (polys->GetNextCell(n, ids))
  {
    if (n < 2) continue;
    for (int i = 0; i < n; ++i)
    {
      vtkIdType a = ids[i];
      vtkIdType b = ids[(i+1)%n];
      uint64_t key = edgeKey(a,b);
      if (edgeMap.find(key) == edgeMap.end())
      {
        double pa[3], pb[3];
        pts->GetPoint(a, pa);
        pts->GetPoint(b, pb);
        double dx = pa[0]-pb[0]; double dy = pa[1]-pb[1]; double dz = pa[2]-pb[2];
        double len = std::sqrt(dx*dx + dy*dy + dz*dz);
        edgeMap[key] = static_cast<vtkIdType>(edges.size());
        edges.push_back({a,b,len});
      }
    }
  }

  // Determine target edge length if not specified: mean edge length
  if (target <= 0.0)
  {
    double sum = 0.0;
    for (const auto& e : edges) sum += e.len;
    if (!edges.empty()) target = sum / static_cast<double>(edges.size());
    else target = 0.0;
  }

  if (target <= 0.0)
  {
    output->ShallowCopy(working);
    return 1;
  }

  // Build point -> neighbor list
  std::vector<std::vector<vtkIdType>> neighbors(np);
  for (const auto& e : edges)
  {
    neighbors[e.a].push_back(e.b);
    neighbors[e.b].push_back(e.a);
  }

  // Locator for projection back to original input surface
  vtkNew<vtkCellLocator> locator;
  locator->SetDataSet(input);
  locator->BuildLocator();

  // Iterative remeshing: split long edges, collapse short edges, then tangential smoothing
  for (int iter = 0; iter < this->MaxIterations; ++iter)
  {
    // Recompute edge lengths
    for (auto &e : edges)
    {
      double pa[3], pb[3];
      pts->GetPoint(e.a, pa);
      pts->GetPoint(e.b, pb);
      double dx = pa[0]-pb[0]; double dy = pa[1]-pb[1]; double dz = pa[2]-pb[2];
      e.len = std::sqrt(dx*dx + dy*dy + dz*dz);
    }

    // Split edges that are too long (> 4/3 * target)
    double splitThresh = (4.0/3.0) * target;
    std::vector<std::tuple<vtkIdType,vtkIdType, vtkIdType>> splits; // a,b,newId
    for (const auto &e : edges)
    {
      if (e.len > splitThresh)
      {
        double pa[3], pb[3];
        pts->GetPoint(e.a, pa);
        pts->GetPoint(e.b, pb);
        double mid[3] = {(pa[0]+pb[0])*0.5, (pa[1]+pb[1])*0.5, (pa[2]+pb[2])*0.5};
        vtkIdType newId = pts->InsertNextPoint(mid);
        splits.emplace_back(e.a, e.b, newId);
      }
    }

    if (!splits.empty())
    {
      // After splits, rebuild neighbors and edges conservatively
      for (auto &s : splits)
      {
        vtkIdType a,b,nid;
        std::tie(a,b,nid) = s;
        neighbors[nid].push_back(a);
        neighbors[nid].push_back(b);
        neighbors[a].push_back(nid);
        neighbors[b].push_back(nid);
      }
    }

    // Collapse short edges (< 2/3 * target)
    double collapseThresh = (2.0/3.0) * target;
    std::vector<char> removed(pts->GetNumberOfPoints(), 0);
    for (const auto &e : edges)
    {
      if (e.len < collapseThresh)
      {
        // collapse b into a (remove b)
        vtkIdType a = e.a; vtkIdType b = e.b;
        if (removed[a] || removed[b]) continue;
        double pa[3], pb[3];
        pts->GetPoint(a, pa);
        pts->GetPoint(b, pb);
        double newp[3] = {(pa[0]+pb[0])*0.5, (pa[1]+pb[1])*0.5, (pa[2]+pb[2])*0.5};
        pts->SetPoint(a, newp);
        removed[b] = 1;
        // rewire neighbors: move b's neighbors to a
        for (vtkIdType nb : neighbors[b])
        {
          if (nb == a) continue;
          auto &vec = neighbors[nb];
          std::replace(vec.begin(), vec.end(), b, a);
          // add nb to a's neighbor list if missing
          if (std::find(neighbors[a].begin(), neighbors[a].end(), nb) == neighbors[a].end())
            neighbors[a].push_back(nb);
        }
        neighbors[b].clear();
      }
    }

    // Compact points by removing removed vertices (lazy: keep them but don't use)

    // Tangential smoothing: for each vertex, move toward Laplacian centroid and then project onto original surface
    vtkNew<vtkPoints> newPts;
    vtkIdType curNp = pts->GetNumberOfPoints();
    newPts->SetNumberOfPoints(curNp);
    for (vtkIdType pid = 0; pid < curNp; ++pid)
    {
      double p[3]; pts->GetPoint(pid, p);
      if (removed[pid] || neighbors[pid].empty())
      {
        newPts->SetPoint(pid, p);
        continue;
      }
      double centroid[3] = {0.0,0.0,0.0};
      for (vtkIdType nb : neighbors[pid])
      {
        double q[3]; pts->GetPoint(nb,q);
        centroid[0] += q[0]; centroid[1] += q[1]; centroid[2] += q[2];
      }
      centroid[0] /= static_cast<double>(neighbors[pid].size());
      centroid[1] /= static_cast<double>(neighbors[pid].size());
      centroid[2] /= static_cast<double>(neighbors[pid].size());

      // Laplacian move
      double lap[3] = { centroid[0]-p[0], centroid[1]-p[1], centroid[2]-p[2] };
      // Move fraction alpha towards laplacian (alpha small for stability)
      double alpha = 0.5;
      double moved[3] = { p[0] + alpha*lap[0], p[1] + alpha*lap[1], p[2] + alpha*lap[2] };

      // Project moved point back onto original surface to preserve shape
      double closest[3]; vtkIdType cellId; int subId; double dist2;
      locator->FindClosestPoint(moved, closest, cellId, subId, dist2);
      newPts->SetPoint(pid, closest);
    }

    working->SetPoints(newPts);
    pts = working->GetPoints();
    // Rebuild edges from updated points (simple rebuild)
    edges.clear(); edgeMap.clear();
    polys = working->GetPolys(); polys->InitTraversal();
    while (polys->GetNextCell(n, ids))
    {
      if (n < 2) continue;
      for (int i = 0; i < n; ++i)
      {
        vtkIdType a = ids[i]; vtkIdType b = ids[(i+1)%n];
        uint64_t key = edgeKey(a,b);
        if (edgeMap.find(key) == edgeMap.end())
        {
          double pa[3], pb[3]; pts->GetPoint(a, pa); pts->GetPoint(b, pb);
          double dx = pa[0]-pb[0]; double dy = pa[1]-pb[1]; double dz = pa[2]-pb[2];
          double len = std::sqrt(dx*dx + dy*dy + dz*dz);
          edgeMap[key] = static_cast<vtkIdType>(edges.size());
          edges.push_back({a,b,len});
        }
      }
    }

    // Rebuild neighbors
    neighbors.clear(); neighbors.resize(pts->GetNumberOfPoints());
    for (const auto& e : edges)
    {
      neighbors[e.a].push_back(e.b);
      neighbors[e.b].push_back(e.a);
    }
  }

  // Finalize output
  output->ShallowCopy(working);
  return 1;
}
