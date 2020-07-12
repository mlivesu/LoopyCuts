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

#ifndef LOOP_MESHER_H
#define LOOP_MESHER_H


#include "vcg/complex/algorithms/create/platonic.h"
#include "vcg/complex/algorithms/update/color.h"

template <class MeshType>
class LoopMesher
{
    typedef typename MeshType::FaceType FaceType;
    typedef typename MeshType::VertexType VertexType;
    typedef typename MeshType::ScalarType ScalarType;
    typedef typename MeshType::CoordType CoordType;

public:

    //    static void Mesh(const MeshType &BaseMesh,
    //                     const std::vector<std::vector<std::pair<size_t,size_t> > > &FaceEdgeLoops,
    //                     ScalarType Width,MeshType &LoopMesh)
    //    {
    //        LoopMesh.Clear();
    //        std::map<CoordType,std::vector<size_t> > PosLoops;
    //        //for each loop
    //        for (size_t i=0;i<FaceEdgeLoops.size();i++)
    //        {
    //            std::vector<CoordType> LoopPos;
    //            for (size_t j=0;j<FaceEdgeLoops[i].size();j++)
    //            {
    //                size_t IndexF=FaceEdgeLoops[i][j].first;
    //                size_t IndexE=FaceEdgeLoops[i][j].second;
    //                CoordType Pos0=BaseMesh.face[IndexF].cP0(IndexE);
    //                CoordType Pos1=BaseMesh.face[IndexF].cP1(IndexE);
    //                LoopPos.push_back(Pos0);
    //                LoopPos.push_back(Pos1);
    //            }
    //            std::sort(LoopPos.begin(),LoopPos.end());
    //            typename std::vector<CoordType>::iterator it;
    //            it = std::unique (LoopPos.begin(), LoopPos.end());
    //            LoopPos.resize( std::distance(LoopPos.begin(),it) );
    //            for (size_t j=0;j<LoopPos.size();j++)
    //                PosLoops[LoopPos[j]].push_back(i);
    //        }

    //        //then sort
    //        typename std::map<CoordType,std::vector<size_t> >::iterator ItPosL;
    //        for (ItPosL=PosLoops.begin();ItPosL!=PosLoops.end();ItPosL++)
    //            std::sort((*ItPosL).second.begin(),(*ItPosL).second.end());

    //        //then mesh
    //        for (size_t i=0;i<FaceEdgeLoops.size();i++)
    //        {
    //            //compute normal
    //            std::map<CoordType,CoordType> VertEdgeN;
    //            for (size_t j=0;j<FaceEdgeLoops[i].size();j++)
    //            {
    //                size_t IndexF=FaceEdgeLoops[i][j].first;
    //                size_t IndexE=FaceEdgeLoops[i][j].second;
    //                CoordType N0=BaseMesh.face[IndexF].cN();
    //                CoordType N1=BaseMesh.face[IndexF].cFFp(IndexE)->cN();
    //                CoordType Pos0=BaseMesh.face[IndexF].cP0(IndexE);
    //                CoordType Pos1=BaseMesh.face[IndexF].cP1(IndexE);
    //                CoordType AddNorm=N0+N1;
    //                AddNorm.Normalize();
    //                if (VertEdgeN.count(Pos0)==0)
    //                    VertEdgeN[Pos0]=AddNorm;
    //                else
    //                    VertEdgeN[Pos0]+=AddNorm;

    //                if (VertEdgeN.count(Pos1)==0)
    //                    VertEdgeN[Pos1]=AddNorm;
    //                else
    //                    VertEdgeN[Pos1]+=AddNorm;
    //            }

    //            for (size_t j=0;j<FaceEdgeLoops[i].size();j++)
    //            {
    //                size_t IndexF=FaceEdgeLoops[i][j].first;
    //                size_t IndexE=FaceEdgeLoops[i][j].second;
    //                CoordType Pos0=BaseMesh.face[IndexF].cP0(IndexE);
    //                CoordType Pos1=BaseMesh.face[IndexF].cP1(IndexE);
    //                CoordType N0=VertEdgeN[Pos0];
    //                CoordType N1=VertEdgeN[Pos1];
    //                N0.Normalize();
    //                N1.Normalize();

    //                //then find the level
    //                assert(PosLoops.count(Pos0)>0);
    //                int Level0=-1;
    //                for (size_t k=0;k<PosLoops[Pos0].size();k++)
    //                    if (PosLoops[Pos0][k]==i){Level0=k;break;}
    //                assert(Level0>=0);

    //                assert(PosLoops.count(Pos1)>0);
    //                int Level1=-1;
    //                for (size_t k=0;k<PosLoops[Pos1].size();k++)
    //                    if (PosLoops[Pos1][k]==i){Level1=k;break;}
    //                assert(Level1>=0);

    //                //then compute the final position
    //                CoordType Pos0Final=Pos0+N0*Level0*Width*1.1*2;
    //                CoordType Pos1Final=Pos1+N1*Level1*Width*1.1*2;
    //                MeshType CylMesh;
    //                vcg::tri::OrientedCylinder(CylMesh,Pos0Final,Pos1Final,Width,false);
    //                vcg::tri::Append<MeshType,MeshType>::Mesh(LoopMesh,CylMesh);
    //            }
    //        }
    //    }

    static void Mesh(MeshType &BaseMesh,
                     const std::vector<std::vector<std::pair<size_t,size_t> > > &FaceEdgeLoops,
                     ScalarType Width,MeshType &LoopMesh)
    {
        vcg::tri::UpdateNormal<MeshType>::PerFaceNormalized(BaseMesh);
        vcg::tri::UpdateNormal<MeshType>::PerVertexAngleWeighted(BaseMesh);
        //vcg::tri::Smooth<MeshType>::VertexNormalLaplacian(BaseMesh,3);

        LoopMesh.Clear();
        std::map<std::pair<CoordType,CoordType>,size_t > CurrLevel;
        std::vector<size_t> LoopLevel;

        //for each loop
        for (size_t i=0;i<FaceEdgeLoops.size();i++)
        {
            LoopLevel.push_back(0);

            for (size_t j=0;j<FaceEdgeLoops[i].size();j++)
            {
                size_t IndexF=FaceEdgeLoops[i][j].first;
                size_t IndexE=FaceEdgeLoops[i][j].second;
                CoordType Pos0=BaseMesh.face[IndexF].cP0(IndexE);
                CoordType Pos1=BaseMesh.face[IndexF].cP1(IndexE);
                std::pair<CoordType,CoordType> Key(std::min(Pos0,Pos1),std::max(Pos0,Pos1));
                if (CurrLevel.count(Key)==0)
                    CurrLevel[Key]=0;
                else
                    CurrLevel[Key]++;

                LoopLevel.back()=std::max(LoopLevel.back(),CurrLevel[Key]);
            }
            //then reupdate the level for all edges
            for (size_t j=0;j<FaceEdgeLoops[i].size();j++)
            {
                size_t IndexF=FaceEdgeLoops[i][j].first;
                size_t IndexE=FaceEdgeLoops[i][j].second;
                CoordType Pos0=BaseMesh.face[IndexF].cP0(IndexE);
                CoordType Pos1=BaseMesh.face[IndexF].cP1(IndexE);
                std::pair<CoordType,CoordType> Key(std::min(Pos0,Pos1),std::max(Pos0,Pos1));
                CurrLevel[Key]=LoopLevel.back();
            }
        }

        MeshType EMesh;

        //then mesh
        for (size_t i=0;i<FaceEdgeLoops.size();i++)
        {
            //            //compute normal
            vcg::Color4b CurrC=vcg::Color4b::Scatter(FaceEdgeLoops.size(),i);
            //                        std::map<CoordType,CoordType> VertEdgeN;
            //                        for (size_t j=0;j<FaceEdgeLoops[i].size();j++)
            //                        {
            //                            size_t IndexF=FaceEdgeLoops[i][j].first;
            //                            size_t IndexE=FaceEdgeLoops[i][j].second;
            //                            CoordType N0=BaseMesh.face[IndexF].cN();
            //                            CoordType N1=BaseMesh.face[IndexF].cFFp(IndexE)->cN();
            //                            CoordType Pos0=BaseMesh.face[IndexF].cP0(IndexE);
            //                            CoordType Pos1=BaseMesh.face[IndexF].cP1(IndexE);
            //                            CoordType AddNorm=N0+N1;
            //                            AddNorm.Normalize();
            //                            if (VertEdgeN.count(Pos0)==0)
            //                                VertEdgeN[Pos0]=AddNorm;
            //                            else
            //                                VertEdgeN[Pos0]+=AddNorm;

            //                            if (VertEdgeN.count(Pos1)==0)
            //                                VertEdgeN[Pos1]=AddNorm;
            //                            else
            //                                VertEdgeN[Pos1]+=AddNorm;
            //                        }

            //find the biggest level
            int CurrLev=LoopLevel[i];
            for (size_t j=0;j<FaceEdgeLoops[i].size();j++)
            {
                size_t IndexF=FaceEdgeLoops[i][j].first;
                size_t IndexE=FaceEdgeLoops[i][j].second;
                CoordType Pos0=BaseMesh.face[IndexF].cP0(IndexE);
                CoordType Pos1=BaseMesh.face[IndexF].cP1(IndexE);
                CoordType N0=BaseMesh.face[IndexF].cV0(IndexE)->cN();//
                CoordType N1=BaseMesh.face[IndexF].cV1(IndexE)->cN();//
                N0.Normalize();
                N1.Normalize();
                //then compute the final position
                CoordType OffSet0=N0*Width*1.1;
                CoordType OffSet1=N1*Width*1.1;
                CoordType Pos0Final=OffSet0+Pos0+N0*CurrLev*Width*1.1*2;
                CoordType Pos1Final=OffSet1+Pos1+N1*CurrLev*Width*1.1*2;
                vcg::tri::Allocator<MeshType>::AddEdge(EMesh,Pos0Final,Pos1Final);
                //EMesh.edge.back().C()=CurrC;
                MeshType CylMesh,Sph0Mesh,Sph1Mesh;
                vcg::tri::OrientedCylinder(CylMesh,Pos0Final,Pos1Final,Width,false);
                vcg::tri::UpdateColor<MeshType>::PerFaceConstant(CylMesh,CurrC);
                vcg::tri::Append<MeshType,MeshType>::Mesh(LoopMesh,CylMesh);

                vcg::tri::Sphere(Sph0Mesh,1);
                for (size_t j=0;j<Sph0Mesh.vert.size();j++)
                {
                    Sph0Mesh.vert[j].P()*=Width;
                    Sph0Mesh.vert[j].P()+=Pos0Final;
                }
                vcg::tri::Sphere(Sph1Mesh,1);
                for (size_t j=0;j<Sph1Mesh.vert.size();j++)
                {
                    Sph1Mesh.vert[j].P()*=Width;
                    Sph1Mesh.vert[j].P()+=Pos1Final;
                }
                vcg::tri::UpdateColor<MeshType>::PerFaceConstant(Sph0Mesh,CurrC);
                vcg::tri::Append<MeshType,MeshType>::Mesh(LoopMesh,Sph0Mesh);
                vcg::tri::UpdateColor<MeshType>::PerFaceConstant(Sph1Mesh,CurrC);
                vcg::tri::Append<MeshType,MeshType>::Mesh(LoopMesh,Sph1Mesh);

            }
            //vcg::tri::Clean<MeshType>::RemoveDuplicateVertex(EMesh);
        }


        //            for (size_t s=0;s<3;s++)
        //            {
        //                std::vector<CoordType> posVec(EMesh.vn,CoordType(0,0,0));
        //                std::vector<int>     cntVec(EMesh.vn,0);

        //                for(int i =0; i<EMesh.en;++i)
        //                {
        //                    for(int j=0;j<2;++j)
        //                    {
        //                        int vertInd = tri::Index(EMesh,EMesh.edge[i].V0(j));
        //                        posVec[vertInd] += EMesh.edge[i].V1(j)->P();
        //                        cntVec[vertInd] += 1;
        //                    }
        //                }
        //                for(int i =0; i<EMesh.vn;++i)
        //                    EMesh.vert[i].P()=posVec[i] /cntVec[i];
        //            }

        //            for(int i =0; i<EMesh.en;++i)
        //            {
        //               MeshType CylMesh,Sph0Mesh,Sph1Mesh;
        //               CoordType Pos0= EMesh.edge[i].P(0);
        //               CoordType Pos1= EMesh.edge[i].P(1);
        ////               vcg::tri::OrientedCylinder(CylMesh, EMesh.edge[i].P(0),EMesh.edge[i].P(1),Width,false);
        ////               vcg::tri::UpdateColor<MeshType>::PerFaceConstant(CylMesh,EMesh.edge[i].C());
        ////               vcg::tri::Append<MeshType,MeshType>::Mesh(LoopMesh,CylMesh);

        //               vcg::tri::Sphere(Sph0Mesh);
        //               for (size_t j=0;j<Sph0Mesh.vert.size();j++)
        //                   Sph0Mesh.vert[j].P()+=Pos0;

        //               vcg::tri::UpdateColor<MeshType>::PerFaceConstant(Sph0Mesh,EMesh.edge[i].C());
        //               vcg::tri::Append<MeshType,MeshType>::Mesh(LoopMesh,Sph0Mesh);

        ////               vcg::tri::Sphere(Sph1Mesh);
        ////               for (size_t j=0;j<Sph1Mesh.vert.size();j++)
        ////                   Sph1Mesh.vert[j].P()+=Pos1;
        ////               vcg::tri::UpdateColor<MeshType>::PerFaceConstant(Sph1Mesh,EMesh.edge[i].C());
        ////               vcg::tri::Append<MeshType,MeshType>::Mesh(LoopMesh,Sph1Mesh);
        //            }


    }
};

#endif
