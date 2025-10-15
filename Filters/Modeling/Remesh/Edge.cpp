#include "Edge.h"
#include "HalfEdge.h"
#include "Vertex.h"

namespace vtkBotschKobbeltRemeshing
{

Edge::Edge()
    : he()
    , index(-1)
    , remove(false)
    , feature(false)
{
}

double Edge::lengthSquared() const
{
    Eigen::Vector3d a = he->vertex->position;
    Eigen::Vector3d b = he->flip->vertex->position;

    return (b - a).squaredNorm();
}

} // namespace vtkBotschKobbeltRemeshing