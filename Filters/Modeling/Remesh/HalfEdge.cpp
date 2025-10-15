#include "HalfEdge.h"

namespace vtkBotschKobbeltRemeshing
{

HalfEdge::HalfEdge()
    : next()
    , flip()
    , vertex()
    , edge()
    , face()
    , onBoundary(false)
    , index(-1)
    , remove(false)
{
}

} // namespace vtkBotschKobbeltRemeshing
