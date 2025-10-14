#ifndef VTKBOTSCHKOBBELTREMESHING_EDGE_H
#define VTKBOTSCHKOBBELTREMESHING_EDGE_H

#include "Types.h"

namespace vtkBotschKobbeltRemeshing
{

class Edge
{
public:
    // one of the two half edges associated with this edge
    HalfEdgeIter he;

    double lengthSquared() const;

    // id between 0 and |E|-1
    int index;

    // flag for deletion
    bool remove;

    // feature flag
    bool feature;
};

} // namespace vtkBotschKobbeltRemeshing

#endif
