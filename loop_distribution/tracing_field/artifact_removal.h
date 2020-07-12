/***************************************************************************/
/* Copyright(C) 2020

 Marco Livesu
 Italian National Research Council

 and

 Nico Pietroni
 University Of Technology Sydney

 All rights reserved.
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
****************************************************************************/

#ifndef ARTIFACT_REMOVAL
#define ARTIFACT_REMOVAL

#include <vcg/simplex/face/topology.h>
#include <vcg/complex/algorithms/update/flag.h>
#include <vcg/space/color4.h>
#include <vcg/complex/algorithms/smooth.h>
#include <vcg/complex/algorithms/parametrization/tangent_field_operators.h>

/**
 * @brief The class that removes small degenerate artifacts from the mesh
 */
template < class MeshType >
class ArtifactRemoval
{
    typedef typename MeshType::CoordType    CoordType;
    typedef typename MeshType::VertexType   VertexType;
    typedef typename MeshType::FaceType     FaceType;
    typedef typename MeshType::ScalarType   ScalarType;

    public:

    static bool IsBadVertex( VertexType &v_test,
                             ScalarType maxDiff=-0.6)
    {
        std::vector<FaceType*> faces;
        std::vector<int> indexes;

        vcg::face::VFStarVF(&v_test,faces,indexes);
        //then check each pair of faces
        for (size_t i=0;i<faces.size()-1;i++)
            for (size_t j=(i+1);j<faces.size();j++)
            {
                CoordType N0=faces[i]->N();
                CoordType N1=faces[j]->N();
                N0.Normalize();
                N1.Normalize();
                if ((N0*N1)<maxDiff)return true;
            }
        return false;
    }

    static void SelectBadVertices(MeshType & mesh,
                                  bool colorize=true)
    {
        vcg::tri::UpdateFlags<MeshType>::VertexClearS(mesh);
        for (size_t i=0;i<mesh.vert.size();i++)
        {
            if (!IsBadVertex(mesh.vert[i]))continue;
            mesh.vert[i].SetS();

            if (!colorize)continue;

            std::vector<FaceType*> faces;
            std::vector<int> indexes;
            vcg::face::VFStarVF(&mesh.vert[i],faces,indexes);
            for (size_t j=0;j<faces.size();j++)
                 faces[j]->C()=vcg::Color4b(255,0,0,255);
        }
    }

    static void SolveBadVertices(MeshType & mesh)
    {
        //first select the bad ones
        SelectBadVertices(mesh);
        //then smooth the selected ones
        vcg::tri::Smooth<MeshType>::VertexCoordLaplacian(mesh,3,true);
        //update per face normals
        //vcg::tri::UpdateNormal<MeshType>::PerVertexPerFace(mesh);
        vcg::tri::UpdateNormal<MeshType>::PerFaceNormalized(mesh);
        vcg::tri::UpdateNormal<MeshType>::PerVertexNormalized(mesh);
        //finally reupdate the tangent field in a coherent way
        vcg::tri::CrossField<MeshType>::AdjustDirectionsOnTangentspace(mesh);

    }
};
#endif //SEPARATRIX_OPTIMIZER
