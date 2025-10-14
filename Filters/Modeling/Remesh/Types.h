#ifndef VTKBOTSCHKOBBELTREMESHING_TYPES_H
#define VTKBOTSCHKOBBELTREMESHING_TYPES_H

#include <cstdlib>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
#include <cmath>

#include "vtk_eigen.h"
#include VTK_EIGEN(Core)
#include VTK_EIGEN(Dense)

namespace vtkBotschKobbeltRemeshing
{

class Vertex;
class Edge;
class Face;
class HalfEdge;
class Mesh;
class MeshIO;

using HalfEdgeIter = std::vector<HalfEdge>::iterator;
using HalfEdgeCIter = std::vector<HalfEdge>::const_iterator;
using VertexIter = std::vector<Vertex>::iterator;
using VertexCIter = std::vector<Vertex>::const_iterator;
using EdgeIter = std::vector<Edge>::iterator;
using EdgeCIter = std::vector<Edge>::const_iterator;
using FaceIter = std::vector<Face>::iterator;
using FaceCIter = std::vector<Face>::const_iterator;
using VectorIter = std::vector<Eigen::Vector3d>::iterator;
using VectorCIter = std::vector<Eigen::Vector3d>::const_iterator;

} // namespace vtkBotschKobbeltRemeshing

#endif
