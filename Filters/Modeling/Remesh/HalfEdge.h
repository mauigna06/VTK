#ifndef VTKBOTSCHKOBBELTREMESHING_HALFEDGE_H
#define VTKBOTSCHKOBBELTREMESHING_HALFEDGE_H

#include "Types.h"

namespace vtkBotschKobbeltRemeshing
{

class HalfEdge
{
public:
    HalfEdge();

    // next halfedge around the current face
    HalfEdgeIter next;

    // other halfedge associated with this edge
    HalfEdgeIter flip;

    // vertex at the tail of the halfedge
    VertexIter vertex;

    // edge associated with this halfedge
    EdgeIter edge;

    // face associated with this halfedge
    FaceIter face;

    // checks if this halfedge is contained in boundary loop
    bool onBoundary;

    // id between 0 and |H|-1
    int index;

    // flag for deletion
    bool remove;
};

} // namespace vtkBotschKobbeltRemeshing

#endif
