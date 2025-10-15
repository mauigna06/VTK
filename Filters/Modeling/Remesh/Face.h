#ifndef VTKBOTSCHKOBBELTREMESHING_FACE_H
#define VTKBOTSCHKOBBELTREMESHING_FACE_H

#include "Types.h"

namespace vtkBotschKobbeltRemeshing
{

class BoundingBox;

class Face
{
public:
    Face();

    // one of the halfedges associated with this face
    HalfEdgeIter he;

    // checks if this face lies on boundary
    bool isBoundary() const;

    // returns face area
    double area() const;

    // returns normal to face
    Eigen::Vector3d normal() const;

    // computes the bounding box of the face
    BoundingBox boundingBox() const;

    // computes the centroid of the face
    Eigen::Vector3d centroid() const;

    // returns closest point on face to p
    double closestPoint(const Eigen::Vector3d& p, Eigen::Vector3d& c) const;

    // id between 0 and |F|-1
    int index;

    // flag for deletion
    bool remove;
};

} // namespace vtkBotschKobbeltRemeshing

#endif
