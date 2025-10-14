#include "Edge.h"
#include "HalfEdge.h"
#include "Vertex.h"

namespace vtkBotschKobbeltRemeshing
{

double Edge::lengthSquared() const
{
    Eigen::Vector3d a = he->vertex->position;
    Eigen::Vector3d b = he->flip->vertex->position;

    return (b - a).squaredNorm();
}

} // namespace vtkBotschKobbeltRemeshing