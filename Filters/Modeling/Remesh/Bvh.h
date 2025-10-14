#ifndef VTKBOTSCHKOBBELTREMESHING_BVH_H
#define VTKBOTSCHKOBBELTREMESHING_BVH_H

#include "Types.h"
#include "BoundingBox.h"

namespace vtkBotschKobbeltRemeshing
{

class Node
{
public:
    // member variables
    BoundingBox boundingBox;
    int startId;
    int range;
    int rightOffset;
};

class Bvh
{
public:
    explicit Bvh(const int leafSize0 = 1);

    // builds the bvh
    void build(Mesh* meshPtr0, BoundingBox& boundingBox);

    // checks if a face overlaps with another face. Returns face id
    double nearestPoint(const double dMin, const Eigen::Vector3d& p, Eigen::Vector3d& np) const;

private:
    int nodeCount;
    int leafCount;
    int leafSize;
    std::vector<Node> flatTree;
    Mesh* meshPtr;
};

} // namespace vtkBotschKobbeltRemeshing

#endif
