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

#ifndef SMOOTH_LOOP_MESH_H
#define SMOOTH_LOOP_MESH_H

#include <vcg/complex/algorithms/update/flag.h>
#include <vcg/complex/algorithms/update/topology.h>
#include <vcg/complex/algorithms/clean.h>
#include <vcg/simplex/face/pos.h>
#include <vcg/space/index/grid_static_obj.h>
#include <vcg/complex/algorithms/closest.h>

class SmoothLoops
{
    typedef typename CMesh::FaceType FaceType;
    typedef typename CMesh::VertexType VertexType;
    typedef typename CMesh::CoordType CoordType;
    typedef typename CMesh::ScalarType ScalarType;


    static void LaplacianStep(CMesh &mesh,
                              std::vector<bool> &MovingV,
                              ScalarType Damp,
                              std::set<std::pair<size_t,size_t> > *ContributeEdge=NULL)
    {
        assert(Damp>=0);
        assert(Damp<=1);
        //cumulate
        std::vector<CoordType> AvgPos(mesh.vert.size(),CoordType(0,0,0));
        std::vector<size_t> NumPos(mesh.vert.size(),0);

        for (size_t i=0;i<mesh.face.size();i++)
        {
            std::vector<size_t> IndexVF;
            for (size_t j=0;j<mesh.face[i].VN();j++)
            {
                size_t IndexV=vcg::tri::Index(mesh,mesh.face[i].V(j));
                IndexVF.push_back(IndexV);
            }

            for (size_t j=0;j<IndexVF.size();j++)
            {
                size_t IndexV0=vcg::tri::Index(mesh,mesh.face[i].V(j));

                CoordType currPos=mesh.vert[IndexV0].P();
                for (size_t j=0;j<IndexVF.size();j++)
                {
                    size_t IndexV1=IndexVF[j];
                    if (IndexV0==IndexV1)continue;//same one

                    std::pair<size_t,size_t> keyE(std::min(IndexV0,IndexV1),
                                                  std::max(IndexV0,IndexV1));
                    if ((ContributeEdge!=NULL)&&
                        ((*ContributeEdge).count(keyE)==0))continue;

                    AvgPos[IndexV1]+=currPos;
                    NumPos[IndexV1]++;
                }
            }
        }
        //average
        for (size_t i=0;i<AvgPos.size();i++)
        {
            if (!MovingV[i])continue;
            if (NumPos[i]==0)continue;
            mesh.vert[i].P()=AvgPos[i]/NumPos[i]*(1-Damp)+mesh.vert[i].P()*Damp;
        }
    }

    static void ReprojectStep(CMesh &mesh,CMesh &target,
                              vcg::GridStaticPtr<FaceType,ScalarType> &Grid)
    {
        for (size_t i=0;i<mesh.vert.size();i++)
        {
            CoordType ClosePt;
            ScalarType maxD,minD;
            maxD=mesh.bbox.Diag();
            FaceType *f=NULL;
            f=vcg::tri::GetClosestFaceBase(target,Grid,mesh.vert[i].P(),maxD,minD,ClosePt);
            assert(f!=NULL);
            mesh.vert[i].P()=ClosePt;
        }
    }

    static void GetFixedVert(CMesh &mesh,
                             std::vector<std::vector<std::pair<size_t,size_t> > > &FaceEdgeLoops,
                             std::vector<bool> &IsVertex,
                             std::vector<bool> &IsLine,
                             std::set<std::pair<size_t,size_t> > &EdgeLoopsSet)
    {
        EdgeLoopsSet.clear();
        std::vector<size_t> NumL(mesh.vert.size(),0);
        for (size_t i=0;i<FaceEdgeLoops.size();i++)
        {
            for (size_t j=0;j<FaceEdgeLoops[i].size();j++)
            {
            size_t IndexF=FaceEdgeLoops[i][j].first;
            size_t IndexE=FaceEdgeLoops[i][j].second;
            size_t indexV0=vcg::tri::Index(mesh,mesh.face[IndexF].V0(IndexE));
            size_t indexV1=vcg::tri::Index(mesh,mesh.face[IndexF].V1(IndexE));
            EdgeLoopsSet.insert(std::pair<size_t,size_t>(std::min(indexV0,indexV1),
                                                         std::max(indexV0,indexV1)));
            NumL[indexV0]++;
            NumL[indexV1]++;
            }
        }

        IsVertex= std::vector<bool>(mesh.vert.size(),false);
        IsLine= std::vector<bool>(mesh.vert.size(),false);
        for (size_t i=0;i<NumL.size();i++)
        {
            if(NumL[i]==0)continue;
            if(NumL[i]==2)IsLine[i]=true;
            else IsVertex[i]=true;
        }
    }

public:

    static void SmoothLaplacian(CMesh &mesh,
                                std::vector<std::vector<std::pair<size_t,size_t> > > FaceEdgeLoops,
                                ScalarType Damp=0.5,
                                int Steps=3)
    {
        vcg::tri::UpdateTopology<CMesh>::FaceFace(mesh);

        vcg::GridStaticPtr<FaceType,ScalarType> Grid;

        CMesh TargetM;
        vcg::tri::Append<CMesh,CMesh>::Mesh(TargetM,mesh);
        Grid.Set(TargetM.face.begin(),TargetM.face.end());


        std::vector<bool> IsVertex,IsLine;
        std::set<std::pair<size_t,size_t> > EdgeLoopsSet;
        GetFixedVert(mesh,FaceEdgeLoops,IsVertex,IsLine,EdgeLoopsSet);

        std::vector<bool> MovingV0=IsLine;
        //std::vector<bool> ContributeV0=IsLine;
        for (size_t i=0;i<MovingV0.size();i++)
        {
            if (IsVertex[i])
            MovingV0[i]=false;
        }

        std::vector<bool> MovingV1=std::vector<bool>(mesh.vert.size(),true);
        //std::vector<bool> ContributeV1=std::vector<bool>(mesh.vert.size(),true);
        for (size_t i=0;i<MovingV1.size();i++)
        {
            if (IsVertex[i])MovingV1[i]=false;
            if (IsLine[i])MovingV1[i]=false;
        }

        for (size_t i=0;i<Steps;i++)
        {
            LaplacianStep(mesh,MovingV0,Damp,&EdgeLoopsSet);
            LaplacianStep(mesh,MovingV1,Damp);
            ReprojectStep(mesh,TargetM,Grid);
        }
    }

};

#endif
