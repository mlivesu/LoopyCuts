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


#ifndef MESH_TYPE_H
#define MESH_TYPE_H

/// vcg imports
#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/complex/algorithms/update/normal.h>
#include <vcg/complex/algorithms/create/platonic.h>
#include <vcg/complex/algorithms/refine.h>
/// wrapper imports
#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/import_field.h>
#include <wrap/gl/trimesh.h>
#include <tracing_field/loop_common_functions.h>
#include <tracing_field/sharp_feature_manager.h>

#include <vcg/complex/algorithms/polygonal_algorithms.h>
#include <vcg/complex/algorithms/parametrization/tangent_field_operators.h>
//#include <vcg/complex/algorithms/curve_on_manifold.h>

using namespace vcg;
class CFace;
class CVertex;
class CEdge;

struct MyLoopUsedTypes : public UsedTypes<	Use<CVertex>::AsVertexType,
        Use<CEdge>::AsEdgeType,
        Use<CFace>::AsFaceType>{};

//compositing wanted proprieties
class CVertex : public vcg::Vertex< MyLoopUsedTypes,
        vcg::vertex::TexCoord2d,
        vcg::vertex::Coord3d,
        vcg::vertex::Normal3d,
        vcg::vertex::BitFlags,
        vcg::vertex::VFAdj,
        vcg::vertex::VEAdj,
        vcg::vertex::Qualityd,
        vcg::vertex::Color4b,
        vcg::vertex::Mark,
        vcg::vertex::CurvatureDird>{};

class CFace   : public vcg::Face<  MyLoopUsedTypes,
        vcg::face::VertexRef,
        vcg::face::Normal3d,
        vcg::face::BitFlags,
        vcg::face::CurvatureDird,
        vcg::face::FFAdj,
        vcg::face::VFAdj,
        vcg::face::Qualityd,
        vcg::face::Color4b,
        vcg::face::Mark,
        vcg::face::CurvatureDird,
        vcg::face::WedgeTexCoord2d>
{
public:
    int IndexFeature;
    int FeaturePath;
    bool Concave;
    int IndexE;

    void ImportData(const CFace  & left )
    {
        vcg::Face<  MyLoopUsedTypes,
                vcg::face::VertexRef,
                vcg::face::Normal3d,
                vcg::face::BitFlags,
                vcg::face::CurvatureDird,
                vcg::face::FFAdj,
                vcg::face::VFAdj,
                vcg::face::Qualityd,
                vcg::face::Color4b,
                vcg::face::Mark,
                vcg::face::CurvatureDird,
                vcg::face::WedgeTexCoord2d>::ImportData(left);

        IndexFeature=left.IndexFeature;
        FeaturePath=left.FeaturePath;
        Concave=left.Concave;
        IndexE=left.IndexE;
    }
};

class CEdge   : public vcg::Edge<  MyLoopUsedTypes,
        vcg::edge::VertexRef,
        vcg::edge::EEAdj,
        vcg::edge::VEAdj,
        vcg::edge::Qualityd,
        vcg::edge::Color4b,
        vcg::edge::VEAdj,
        vcg::edge::BitFlags> {};



class CMesh   : public vcg::tri::TriMesh< std::vector<CVertex>,std::vector<CEdge>,std::vector<CFace> >
{
public:

    //ScalarType FlatDegree;

    void UpdateAttributes()
    {
        vcg::tri::UpdateNormal<CMesh>::PerFaceNormalized(*this);
        vcg::tri::UpdateNormal<CMesh>::PerVertexNormalized(*this);
        vcg::tri::UpdateBounding<CMesh>::Box(*this);
        vcg::tri::UpdateTopology<CMesh>::FaceFace(*this);
        vcg::tri::UpdateTopology<CMesh>::VertexFace(*this);
        vcg::tri::UpdateFlags<CMesh>::FaceBorderFromFF(*this);
        vcg::tri::UpdateFlags<CMesh>::VertexBorderFromFaceBorder(*this);
    }

    bool LoadField(std::string field_filename)
    {
        int position0=field_filename.find(".ffield");
        int position1=field_filename.find(".rosy");


        if (position0!=-1)
        {
            bool loaded=vcg::tri::io::ImporterFIELD<CMesh>::LoadFFIELD(*this,field_filename.c_str());
            if (!loaded)return false;
            vcg::tri::CrossField<CMesh>::OrientDirectionFaceCoherently(*this);
            vcg::tri::CrossField<CMesh>::UpdateSingularByCross(*this);
            return true;
        }
        if (position1!=-1)
        {
            std::cout<<"Importing ROSY field"<<std::endl;
            bool loaded=vcg::tri::io::ImporterFIELD<CMesh>::Load4ROSY(*this,field_filename.c_str());
            std::cout<<"Imported ROSY field"<<std::endl;
            if (!loaded)return false;
            vcg::tri::CrossField<CMesh>::OrientDirectionFaceCoherently(*this);
            vcg::tri::CrossField<CMesh>::UpdateSingularByCross(*this);
            return true;
        }
        return false;
    }

    bool LoadFromFile(std::string filename)
    {
        Clear();
        if(filename.empty()) return false;
        int position0=filename.find(".ply");
        int position1=filename.find(".obj");
        int position2=filename.find(".off");

        if (position0!=-1)
        {
            int err=vcg::tri::io::ImporterPLY<CMesh>::Open(*this,filename.c_str());
            if (err!=vcg::ply::E_NOERROR)return false;
            return true;
        }
        if (position1!=-1)
        {
            int mask;
            vcg::tri::io::ImporterOBJ<CMesh>::LoadMask(filename.c_str(),mask);
            int err=vcg::tri::io::ImporterOBJ<CMesh>::Open(*this,filename.c_str(),mask);
            if ((err!=0)&&(err!=5))return false;
            return true;
        }
        if (position2!=-1)
        {
            int err=vcg::tri::io::ImporterOFF<CMesh>::Open(*this,filename.c_str());
            if (err!=0)return false;
            return true;
        }
        return false;
    }


    size_t NumEdgeSel()
    {
        size_t NumCrease=0;
        for (size_t i=0;i<face.size();i++)
            for (size_t j=0;j<3;j++)
                if (face[i].IsFaceEdgeS(j))NumCrease++;
        return NumCrease;
    }


    void GLDrawSharpEdges()
    {
        glPushAttrib(GL_ALL_ATTRIB_BITS);
        glDisable(GL_LIGHTING);
        glDepthRange(0,0.9999);
        glLineWidth(10);
        glBegin(GL_LINES);
        for (size_t i=0;i<face.size();i++)
            for (size_t j=0;j<3;j++)
            {
                if (!face[i].IsFaceEdgeS(j))continue;

                vcg::glColor(vcg::Color4b(255,0,0,255));

                CoordType Pos0=face[i].P0(j);
                CoordType Pos1=face[i].P1(j);
                vcg::glVertex(Pos0);
                vcg::glVertex(Pos1);
            }
        glEnd();
        glPopAttrib();
    }
};

#endif
