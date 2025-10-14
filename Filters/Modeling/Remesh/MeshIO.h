#ifndef VTKBOTSCHKOBBELTREMESHING_MESHIO_H
#define VTKBOTSCHKOBBELTREMESHING_MESHIO_H

#include <fstream>

#include "Types.h"

class vtkPolyData;

namespace vtkBotschKobbeltRemeshing
{

class MeshData;

class MeshIO
{
public:
    // reads data from obj file
    static bool read(std::ifstream& in, Mesh& mesh);

    // writes data in obj format
    static void write(std::ofstream& out, const Mesh& mesh);

    // sets index for vertices
    static void indexElements(Mesh& mesh);

    // builds mesh from vtkPolyData
    static bool BuildFromPolyData(vtkPolyData* poly, Mesh& mesh);

    // copies mesh content into vtkPolyData
    static void CopyToPolyData(const Mesh& mesh, vtkPolyData* poly);

private:
    // reserves spave for mesh vertices, uvs, normals and faces
    static void preallocateMeshElements(const MeshData& data, Mesh& mesh);

    // checks if any vertex is not contained in a face
    static void checkIsolatedVertices(const Mesh& mesh);

    // checks if a vertex is non-manifold
    static void checkNonManifoldVertices(const Mesh& mesh);

    // builds the halfedge mesh
    static bool buildMesh(const MeshData& data, Mesh& mesh);
};

} // namespace vtkBotschKobbeltRemeshing

#endif
