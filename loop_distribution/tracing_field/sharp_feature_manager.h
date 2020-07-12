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


#ifndef SHARP_FEATURE_MANAGER
#define SHARP_FEATURE_MANAGER

//vcg imports
#include <vcg/complex/complex.h>
//#include <vcg/complex/algorithms/update/bounding.h>
//#include <vcg/complex/algorithms/update/normal.h>
//#include <vcg/complex/algorithms/create/platonic.h>
//#include <vcg/complex/algorithms/refine.h>
//wrapper imports
//#include <wrap/io_trimesh/import.h>
//#include <wrap/io_trimesh/import_field.h>
#ifndef NO_TRACING_OPENGL
#include <wrap/gl/trimesh.h>
#endif
#include <tracing_field/loop_common_functions.h>
//#include <vcg/complex/algorithms/polygonal_algorithms.h>
//#include <vcg/complex/algorithms/parametrization/tangent_field_operators.h>
//#include <vcg/complex/algorithms/curve_on_manifold.h>
//#include <wrap/gl/trimesh.h>

template <class MeshType>
size_t MarkEdgeConnectedComponents(MeshType &mesh,bool AvoidCrossS=false)
{
    //set the quality as -1
    for (size_t i=0;i<mesh.edge.size();i++)
        mesh.edge[i].Q()=-1;

    vcg::tri::UpdateTopology<MeshType>::EdgeEdge(mesh);
    size_t currF=0;
    for (size_t i=0;i<mesh.edge.size();i++)
    {
        if (mesh.edge[i].Q()!=-1)continue;

        std::vector<int> EdgeStack;
        EdgeStack.push_back(i);
        do{
            int CurrE=EdgeStack.back();
            EdgeStack.pop_back();
            if (mesh.edge[CurrE].Q()!=-1)continue;
            mesh.edge[CurrE].Q()=currF;

            int NextI0=vcg::tri::Index(mesh,mesh.edge[CurrE].EEp(0));
            if ((AvoidCrossS)&&(mesh.edge[CurrE].V(0)->IsS()))
                NextI0=-1;
            int NextI1=vcg::tri::Index(mesh,mesh.edge[CurrE].EEp(1));
            if ((AvoidCrossS)&&(mesh.edge[CurrE].V(1)->IsS()))
                NextI1=-1;

            if ((NextI0!=(int)i)&&(NextI0>=0))EdgeStack.push_back(NextI0);
            if ((NextI1!=(int)i)&&(NextI1>=0))EdgeStack.push_back(NextI1);
        }while (!EdgeStack.empty());
        currF++;
    }
    return currF;
    //std::cout<<"Found "<<currF<<" components"<<std::endl;
}

template <class MeshType>
class SharpFeaturesManager
{
    typedef typename MeshType::FaceType FaceType;
    typedef typename MeshType::VertexType VertexType;
    typedef typename MeshType::CoordType CoordType;
    typedef typename MeshType::ScalarType ScalarType;

public:

    //ScalarType FlatDegree;

    //FUNCTIONS ON SEQUENCES OF EDGES (SHARP FEATURES)
    std::vector<std::vector<vcg::face::Pos<FaceType> > > FeatureSeq;
    std::vector<SharpFeatureType> FeatureType;

    //std::vector<bool> ProblematicSeq;
    std::vector<std::vector<size_t> > IndexF0,IndexF1;
    std::vector<int> PreferredDir;

    //std::set<std::pair<CoordType,CoordType> > SharpFeaturesGeo;

    MeshType &mesh;
    SharpFeaturesManager(MeshType &_mesh):mesh(_mesh){}


    typedef std::pair<CoordType,CoordType> EdgeCoordKey;

    std::map<EdgeCoordKey,SharpFeatureType> PerSharpType;

    //    //EDGE SPLITTING FUNCTIONS TO COMPLAIN WITH MULTIPLE FEATURES
    //    void InitFeatureSeqFromEdgeSel(size_t MaxComp,MeshType &Feature)
    //    {
    //        FeatureSeq.clear();
    //        FeatureType.clear();

    //        FeatureSeq.resize(MaxComp);

    //        std::map<std::pair<CoordType,CoordType>,std::pair<int,int> > CoordFaceEdge;
    //        std::map<std::pair<CoordType,CoordType>,SharpFeatureType > CoordEdgeType;
    //        for (size_t i=0;i<mesh.face.size();i++)
    //            for (size_t j=0;j<3;j++)
    //            {
    //                if (mesh.face[i].V0(j)>mesh.face[i].V1(j))continue;
    //                //if (face[i].IsF(j))continue;
    //                if (!mesh.face[i].IsFaceEdgeS(j))continue;

    //                CoordType Pos0=mesh.face[i].P0(j);
    //                CoordType Pos1=mesh.face[i].P1(j);

    //                //then add the features
    //                std::pair<CoordType,CoordType> Key(std::min(Pos0,Pos1),std::max(Pos0,Pos1));

    //                CoordFaceEdge[Key]=std::pair<int,int>(i,j);

    ////                if (IsConcaveEdge(mesh.face[i],j))CoordEdgeType[Key]=ConcaveEdge;
    ////                if (!IsConcaveEdge(mesh.face[i],j))CoordEdgeType[Key]=ConvexEdge;
    //                if (IsConcaveEdgeByTable(mesh.face[i],j))CoordEdgeType[Key]=ConcaveEdge;
    //                if (!IsConcaveEdgeByTable(mesh.face[i],j))CoordEdgeType[Key]=ConvexEdge;
    //            }

    //        for (size_t i=0;i<Feature.edge.size();i++)
    //        {
    //            CoordType Pos0=Feature.edge[i].P(0);
    //            CoordType Pos1=Feature.edge[i].P(1);

    //            std::pair<CoordType,CoordType> Key(std::min(Pos0,Pos1),std::max(Pos0,Pos1));
    //            //get the index of the sequence
    //            size_t IndexComp=Feature.edge[i].Q();
    //            assert(IndexComp>=0);
    //            assert(IndexComp<MaxComp);
    //            assert(CoordFaceEdge.count(Key)>0);
    //            int FaceI=CoordFaceEdge[Key].first;
    //            int EdgeI=CoordFaceEdge[Key].second;
    //            FeatureSeq[IndexComp].push_back(vcg::face::Pos<FaceType>(&mesh.face[FaceI],EdgeI));
    //        }

    //        //set the feature type
    //        FeatureType.resize(MaxComp,ConvexEdge);
    //        //ProblematicSeq.resize(FeatureSeq.size(),false);
    //        for (size_t i=0;i<FeatureSeq.size();i++)
    //        {
    //            std::vector<vcg::face::Pos<FaceType> > ConvexE;
    //            std::vector<vcg::face::Pos<FaceType> > ConcaveE;

    //            for (size_t j=0;j<FeatureSeq[i].size();j++)
    //            {
    //                if (IsConcaveEdge(*FeatureSeq[i][j].F(),FeatureSeq[i][j].E()))
    //                    ConcaveE.push_back(FeatureSeq[i][j]);
    //                else
    //                    ConvexE.push_back(FeatureSeq[i][j]);
    //            }
    //            //first case, all convex
    //            if ((ConcaveE.size()==0)&&(ConvexE.size()>0))
    //                FeatureType[i]=ConvexEdge;
    //            //second case, all concave
    //            if ((ConcaveE.size()>0)&&(ConvexE.size()==0))
    //                FeatureType[i]=ConcaveEdge;
    //            //third case ,force to concave
    //            if ((ConcaveE.size()>0)&&(ConvexE.size()>0))
    //            {
    //                //ProblematicSeq[i]=true;
    //                std::cout<<"WARNIGN NON SPLITTED CHANGE OF CONCAVITY"<<std::endl;
    //                assert(0);
    ////                FeatureType[i]=ConvexEdge;
    ////                FeatureSeq[i]=ConvexE;
    ////                FeatureSeq[i].insert(FeatureSeq[i].end(),ConcaveE.begin(),ConcaveE.end());
    //            }
    //        }
    //    }

    //EDGE SPLITTING FUNCTIONS TO COMPLAIN WITH MULTIPLE FEATURES IT GET THE
    //NUMBER OF FEATURES AND THE EDGE MESH
    void InitFeatureSeqFromEdgeSel(size_t MaxComp,MeshType &Feature)
    {
        FeatureSeq.clear();
        FeatureType.clear();

        FeatureSeq.resize(MaxComp);

        std::map<std::pair<CoordType,CoordType>,std::pair<int,int> > CoordFaceEdge;
        std::map<std::pair<CoordType,CoordType>,SharpFeatureType > CoordEdgeType;
        for (size_t i=0;i<mesh.face.size();i++)
            for (size_t j=0;j<3;j++)
            {
                if (mesh.face[i].V0(j)>mesh.face[i].V1(j))continue;
                //if (face[i].IsF(j))continue;
                if (!mesh.face[i].IsFaceEdgeS(j))continue;

                CoordType Pos0=mesh.face[i].P0(j);
                CoordType Pos1=mesh.face[i].P1(j);

                //then add the features
                std::pair<CoordType,CoordType> Key(std::min(Pos0,Pos1),std::max(Pos0,Pos1));

                CoordFaceEdge[Key]=std::pair<int,int>(i,j);

                //                if (IsConcaveEdge(mesh.face[i],j))CoordEdgeType[Key]=ConcaveEdge;
                //                if (!IsConcaveEdge(mesh.face[i],j))CoordEdgeType[Key]=ConvexEdge;
                if (IsConcaveEdgeByTable(mesh.face[i],j))CoordEdgeType[Key]=ConcaveEdge;
                if (!IsConcaveEdgeByTable(mesh.face[i],j))CoordEdgeType[Key]=ConvexEdge;
            }

        for (size_t i=0;i<Feature.edge.size();i++)
        {
            CoordType Pos0=Feature.edge[i].P(0);
            CoordType Pos1=Feature.edge[i].P(1);

            std::pair<CoordType,CoordType> Key(std::min(Pos0,Pos1),std::max(Pos0,Pos1));
            //get the index of the sequence
            size_t IndexComp=Feature.edge[i].Q();
            assert(IndexComp>=0);
            assert(IndexComp<MaxComp);
            assert(CoordFaceEdge.count(Key)>0);
            int FaceI=CoordFaceEdge[Key].first;
            int EdgeI=CoordFaceEdge[Key].second;
            assert(FaceI<(int)mesh.face.size());
            assert(EdgeI<3);
            assert(EdgeI>=0);
            FeatureSeq[IndexComp].push_back(vcg::face::Pos<FaceType>(&mesh.face[FaceI],EdgeI));
        }

        //set the feature type, by default convex
        FeatureType.resize(MaxComp,ConvexEdge);

        for (size_t i=0;i<FeatureSeq.size();i++)
        {
            std::vector<vcg::face::Pos<FaceType> > ConvexE;
            std::vector<vcg::face::Pos<FaceType> > ConcaveE;

            for (size_t j=0;j<FeatureSeq[i].size();j++)
            {
                if (IsConcaveEdgeByTable(*FeatureSeq[i][j].F(),FeatureSeq[i][j].E()))
                    ConcaveE.push_back(FeatureSeq[i][j]);
                else
                    ConvexE.push_back(FeatureSeq[i][j]);
            }
            //first case, all convex
            if ((ConcaveE.size()==0)&&(ConvexE.size()>0))
                FeatureType[i]=ConvexEdge;
            //second case, all concave
            if ((ConcaveE.size()>0)&&(ConvexE.size()==0))
                FeatureType[i]=ConcaveEdge;
            //third case ,force to concave
            if ((ConcaveE.size()>0)&&(ConvexE.size()>0))
            {
                //ProblematicSeq[i]=true;
                std::cout<<"WARNIGN NON SPLITTED CHANGE OF CONCAVITY"<<std::endl;
                assert(0);
                //                FeatureType[i]=ConvexEdge;
                //                FeatureSeq[i]=ConvexE;
                //                FeatureSeq[i].insert(FeatureSeq[i].end(),ConcaveE.begin(),ConcaveE.end());
            }
        }
    }

    //SET THE PREFERRED DIR FOR EACH FEATURE: USED TO TRACE CONCAVE FEATURES
    void InitPreferredDir()
    {
        //std::cout<<"Creases Sel 2 "<<NumEdgeSel()<<std::endl;

        IndexF0.clear();
        IndexF1.clear();
        PreferredDir.clear();
        //std::cout<<"Creases Sel 3 "<<NumEdgeSel()<<std::endl;

        LoopFunctions<MeshType>::GetSharpFacesInfo(mesh,FeatureSeq,IndexF0,IndexF1,PreferredDir);

        //std::cout<<"Creases Sel 4 "<<NumEdgeSel()<<std::endl;

    }

    //SELECT THE VERTIECS THAT NEED TO BE SPLITTED IN THE SHARP FEATURE SEQUENCE
    //CONSIDERIGN MULTIPLE SHARP FEATURES ON IT
    void SelectSplitByMultipleComponents()
    {
        std::vector<size_t> PerVertSharp(mesh.vert.size(),0);

        for (size_t i=0;i<mesh.face.size();i++)
        {
            for (size_t j=0;j<3;j++)
            {
                //if (face[i].IsF(j))continue;
                if (!mesh.face[i].IsFaceEdgeS(j))continue;

                size_t IndexV0=vcg::tri::Index(mesh,mesh.face[i].V0(j));
                size_t IndexV1=vcg::tri::Index(mesh,mesh.face[i].V1(j));

                PerVertSharp[IndexV0]++;
                PerVertSharp[IndexV1]++;
            }
        }

        //then select the one that join multiple sharp features
        for (size_t i=0;i<PerVertSharp.size();i++)
        {
            assert(PerVertSharp[i]%2==0);
            if (PerVertSharp[i]==0)continue;//no sharp features
            if (PerVertSharp[i]==4)continue;//regular vertices
            mesh.vert[i].SetS();
        }
    }

    void InitFaceEdgeSelFromEdgeSeq()
    {
        vcg::tri::UpdateFlags<MeshType>::FaceClearFaceEdgeS(mesh);
        for (size_t i=0;i<FeatureSeq.size();i++)
            for (size_t j=0;j<FeatureSeq[i].size();j++)
            {
                vcg::face::Pos<FaceType> currPos=FeatureSeq[i][j];
                FaceType *currF=currPos.F();
                size_t currE=currPos.E();
                currPos.FlipF();
                FaceType *Fopp=currPos.F();
                size_t oppE=currPos.E();
                currF->SetFaceEdgeS(currE);
                Fopp->SetFaceEdgeS(oppE);
            }
        //std::cout<<"Creases "<<NumCrease<<std::endl;
    }

    //    void InitFaceEdgeSelFromFeatureGeo()
    //    {
    //        size_t NumCrease=0;
    //        vcg::tri::UpdateFlags<MeshType>::FaceClearFaceEdgeS(mesh);
    //        for (size_t i=0;i<mesh.face.size();i++)
    //            for (size_t j=0;j<3;j++)
    //            {
    //                CoordType Pos0=mesh.face[i].P0(j);
    //                CoordType Pos1=mesh.face[i].P1(j);

    //                //then add the features
    //                std::pair<CoordType,CoordType> Key(std::min(Pos0,Pos1),std::max(Pos0,Pos1));
    //                if(SharpFeaturesGeo.count(Key)==0)continue;
    //                mesh.face[i].SetFaceEdgeS(j);

    //                NumCrease++;
    //            }
    //        std::cout<<"Creases "<<NumCrease<<std::endl;
    //    }

    //        bool IsConcaveEdgeByTable(const FaceType &f0,int IndexE)
    //        {
    //            FaceType *f1=f0.cFFp(IndexE);
    //            if (f1==&f0)return false;
    //            CoordType N0=f0.cN();
    //            CoordType N1=f1->cN();
    //            CoordType EdgeDir=f0.cP1(IndexE)-f0.cP0(IndexE);
    //            EdgeDir.Normalize();
    //            CoordType Cross=N0^N1;
    //            return ((Cross*EdgeDir)<0);
    //        }

    bool IsConcaveEdgeByTable(const FaceType &f0,int IndexE)
    {
        CoordType Pos0=f0.cP0(IndexE);
        CoordType Pos1=f0.cP1(IndexE);
        //then add the features
        EdgeCoordKey Key(std::min(Pos0,Pos1),std::max(Pos0,Pos1));
        assert(PerSharpType.count(Key)>0);
        return (PerSharpType[Key]==ConcaveEdge);
    }

    //SELECT THE VERTIECS THAT NEED TO BE SPLITTED IN THE SHARP FEATURE SEQUENCE
    //CONSIDERIGN CHANGE OF CONCAVITY
    void SelectSplitByConcavityChange()
    {
        std::vector<std::vector<SharpFeatureType> > PerFaceSharp;
        std::vector<std::vector<SharpFeatureType> > PerVertSharp;

        PerFaceSharp.resize(mesh.face.size(),std::vector<SharpFeatureType>(3,NoSharp));

        //set per face features
        for (size_t i=0;i<mesh.face.size();i++)
        {
            //check the number of faux
            //int NumF=0;
            for (size_t j=0;j<3;j++)
            {
                //if (face[i].IsF(j))continue;
                if (!mesh.face[i].IsFaceEdgeS(j))continue;

                //                if (IsConcaveEdge(mesh.face[i],j))
                //                    PerFaceSharp[i][j]=ConcaveEdge;
                //                else
                //                    PerFaceSharp[i][j]=ConvexEdge;
                if (IsConcaveEdgeByTable(mesh.face[i],j))
                    PerFaceSharp[i][j]=ConcaveEdge;
                else
                    PerFaceSharp[i][j]=ConvexEdge;
                //NumF++;
            }
            // assert(NumF<=1);
        }

        //then cumulate per vertex
        PerVertSharp.resize(mesh.vert.size());
        for (size_t i=0;i<mesh.face.size();i++)
        {
            for (size_t j=0;j<3;j++)
            {
                //if (face[i].IsF(j))continue;
                if (!mesh.face[i].IsFaceEdgeS(j))continue;

                SharpFeatureType currFeature=PerFaceSharp[i][j];

                if (currFeature==NoSharp)continue;

                size_t IndexV0=vcg::tri::Index(mesh,mesh.face[i].V0(j));
                size_t IndexV1=vcg::tri::Index(mesh,mesh.face[i].V1(j));

                PerVertSharp[IndexV0].push_back(currFeature);
                PerVertSharp[IndexV1].push_back(currFeature);
            }
        }

        //then select the one that have concavity change
        for (size_t i=0;i<PerVertSharp.size();i++)
        {
            if (PerVertSharp[i].size()==0)continue;//no sharp features
            SharpFeatureType FeatureType0=PerVertSharp[i][0];
            for (size_t j=1;j<PerVertSharp[i].size();j++)
            {
                if (PerVertSharp[i][j]!=FeatureType0)
                {
                    mesh.vert[i].SetS();
                    break;
                }
            }
        }
    }

    size_t NumVertSel()
    {
        size_t NumVSel=0;
        for (size_t i=0;i<mesh.vert.size();i++)
            if (mesh.vert[i].IsS())NumVSel++;
        return NumVSel;
    }

    void ColorByPartition()
    {
        vcg::tri::UpdateColor<MeshType>::PerFaceConstant(mesh);

        size_t Range=(IndexF0.size()+IndexF1.size());
        assert(IndexF0.size()==IndexF1.size());
        for (size_t i=0;i<IndexF0.size();i++)
        {
            if (FeatureType[i]!=ConcaveEdge)continue;
            int IndexP0=i*2;
            int IndexP1=i*2+1;
            vcg::Color4b currCol0=vcg::Color4b::Scatter(Range,IndexP0);
            vcg::Color4b currCol1=vcg::Color4b::Scatter(Range,IndexP1);
            for (size_t j=0;j<IndexF0[i].size();j++)
            {
                size_t IndexF=IndexF0[i][j];
                mesh.face[IndexF].C()=currCol0;
            }
            for (size_t j=0;j<IndexF1[i].size();j++)
            {
                size_t IndexF=IndexF1[i][j];
                mesh.face[IndexF].C()=currCol1;
            }
        }
    }

    void SplitByDirectionChange()
    {
        std::cout<<"PLEASE DO SPLIT BY DIRECTION CHANGE"<<std::endl;
    }

    //    void SelectSplitByDirChange(const std::vector<vcg::face::Pos<FaceType> > &ExclDir)
    //    {
    //        typedef std::pair<size_t,size_t> FaceDir;
    //        std::vector<std::vector<FaceDir> > PerFaceDirs;
    //        std::vector<std::vector<FaceDir> > PerVertDirs;

    //        std::set<std::pair<size_t,size_t> > ExclSet;
    //        for (size_t i=0;i<ExclDir.size();i++)
    //        {
    //            vcg::face::Pos<FaceType> CurrPos=ExclDir[i];
    //            size_t IndexV0=vcg::tri::Index(mesh,CurrPos.V());
    //            size_t IndexV1=vcg::tri::Index(mesh,CurrPos.VFlip());
    //            ExclSet.insert(std::pair<size_t,size_t>(std::min(IndexV0,IndexV1),std::max(IndexV0,IndexV1)));
    //        }
    //        //        std::vector<std::vector<SharpFeatureType> > PerFaceSharp;
    //        //        std::vector<std::vector<SharpFeatureType> > PerVertSharp;

    //        //second cycle set per face feature
    //        PerFaceDirs.resize(mesh.face.size());
    //        //        PerFaceSharp.resize(face.size(),std::vector<SharpFeatureType>(3,NoSharp));

    //        for (size_t i=0;i<mesh.face.size();i++)
    //        {
    //            //check the number of faux
    //            for (size_t j=0;j<3;j++)
    //            {
    //                //if (face[i].IsF(j))continue;
    //                if (!mesh.face[i].IsFaceEdgeS(j))continue;
    //                size_t IndexV0=vcg::tri::Index(mesh,mesh.face[i].V0(j));
    //                size_t IndexV1=vcg::tri::Index(mesh,mesh.face[i].V1(j));
    //                std::pair<size_t,size_t> Key(std::min(IndexV0,IndexV1),std::max(IndexV0,IndexV1));
    //                if (ExclSet.count(Key)>0)continue;

    //                vcg::face::Pos<FaceType> CurrPos(&mesh.face[i],j);
    //                int IndexDir=LoopFunctions<MeshType>::getFaceEdgeOrientedDir(CurrPos);
    //                PerFaceDirs[i].push_back(FaceDir(i,IndexDir));

    //                //                if (IsConcaveEdge(face[i],j))
    //                //                    PerFaceSharp[i][j]=ConcaveEdge;
    //                //                else
    //                //                    PerFaceSharp[i][j]=ConvexEdge;

    //                break;
    //            }
    //        }

    //        //thirt cycle, cumulate on vertices
    //        PerVertDirs.resize(mesh.vert.size());
    //        //       PerVertSharp.resize(vert.size());
    //        for (size_t i=0;i<mesh.face.size();i++)
    //        {
    //            for (size_t j=0;j<3;j++)
    //            {
    //                //if (face[i].IsF(j))continue;
    //                if (!mesh.face[i].IsFaceEdgeS(j))continue;

    //                //                vcg::face::Pos<FaceType> Pos0(&face[i],j);
    //                //                vcg::face::Pos<FaceType> Pos1=Pos0;
    //                //                Pos1.FlipF();

    //                //                if (ExclSet.count(Pos0)>0){continue;
    //                //                if (ExclSet.count(Pos1)>0)continue;
    //                //                SharpFeatureType currFeature=PerFaceSharp[i][j];

    //                size_t IndexV0=vcg::tri::Index(mesh,mesh.face[i].V0(j));
    //                size_t IndexV1=vcg::tri::Index(mesh,mesh.face[i].V1(j));

    //                PerVertDirs[IndexV0].insert(PerVertDirs[IndexV0].end(),PerFaceDirs[i].begin(),PerFaceDirs[i].end());
    //                PerVertDirs[IndexV1].insert(PerVertDirs[IndexV1].end(),PerFaceDirs[i].begin(),PerFaceDirs[i].end());

    //                //                if (currFeature==NoSharp)continue;

    //                //                PerVertSharp[IndexV0].push_back(currFeature);
    //                //                PerVertSharp[IndexV1].push_back(currFeature);
    //            }
    //        }

    //        for (size_t i=0;i<PerVertDirs.size();i++)
    //        {
    //            if (mesh.vert[i].IsS())continue;//already selected

    //            if (PerVertDirs[i].size()==0)continue;

    //            //check if change of direction
    //            int IndexF0=PerVertDirs[i][0].first;
    //            int IndexDir0=PerVertDirs[i][0].second;

    //            //            //if no convex no need to split
    //            //            SharpFeatureType FeatureType0=PerVertSharp[i][0];
    //            //            if (FeatureType0==ConvexEdge)continue;

    //            for (size_t j=0;j<PerVertDirs[i].size();j++)
    //            {
    //                int IndexF1=PerVertDirs[i][j].first;
    //                int IndexDir1=PerVertDirs[i][j].second;
    //                int DirTest1=vcg::tri::CrossField<MeshType>::FollowDirection(mesh.face[IndexF0],mesh.face[IndexF1],IndexDir0);
    //                if ((IndexDir1 % 2)!=(DirTest1 % 2))
    //                {
    //                    mesh.vert[i].SetS();
    //                    break;
    //                }
    //            }
    //        }
    //    }


    //    void UpdateEdgeSeqByFieldChange()
    //    {
    //        vcg::tri::UpdateFlags<MeshType>::FaceClearS(mesh);
    //        vcg::tri::UpdateFlags<MeshType>::VertexClearS(mesh);

    //        InitPreferredDir();

    //        //then get the cicles
    //        std::vector<vcg::face::Pos<FaceType> > ExcludeEdges;
    //        for (size_t i=0;i<FeatureSeq.size();i++)
    //        {
    //            if (FeatureType[i]==ConvexEdge)
    //            {
    //                ExcludeEdges.insert(ExcludeEdges.end(),FeatureSeq[i].begin(),FeatureSeq[i].end());
    //            }else
    //                if (LoopFunctions<MeshType>::IsClosedLoop(FeatureSeq[i]))
    //                    ExcludeEdges.insert(ExcludeEdges.end(),FeatureSeq[i].begin(),FeatureSeq[i].end());
    //        }
    //        std::cout<<"exclude "<<ExcludeEdges.size()<<" convex/closed loops edeges from change of directions "<<std::endl;

    //        //std::cout<<"2 Selected "<<NumVertSel()<<" vertices "<<std::endl;
    //        InitFaceEdgeSelFromFeatureGeo();
    //        SelectSplitByMultipleComponents();
    //        SelectSplitByConcavityChange();
    //        SelectSplitByDirChange(ExcludeEdges);

    //        //vcg::tri::io::ExporterPLY<MeshType>::Save(*this,"testSel.ply",vcg::tri::io::Mask::IOM_VERTFLAGS);

    //        std::cout<<"3 Selected "<<NumVertSel()<<" vertices "<<std::endl;

    //        MeshType Feature;
    //        CreateEdgeMeshFromEdgeSel(Feature);
    //        std::cout<<"There are Sharp edge mesh 1:"<<Feature.edge.size()<<std::endl;
    //        //size_t MaxComp=Feature.SplitEdgeConnectedComponents(true);
    //        size_t MaxComp=SplitEdgeConnectedComponents(Feature,true);

    //        InitFeatureSeqFromEdgeSel(MaxComp,Feature);

    //        std::cout<<"There are Sharp sequences:"<<FeatureSeq.size()<<std::endl;

    //        InitPreferredDir();
    //        ColorByPartition();
    //    }

    //    //EDGE SPLITTING FUNCTIONS TO COMPLAIN WITH MULTIPLE FEATURES
    //    typedef std::pair<CoordType,CoordType> EdgeCoordKey;

    //    // Basic subdivision class
    //    struct SplitLev : public   std::unary_function<vcg::face::Pos<FaceType> ,CoordType >
    //    {
    //        std::map<EdgeCoordKey,CoordType> *SplitOps;

    //        void operator()(VertexType &nv,vcg::face::Pos<FaceType>  ep)
    //        {
    //            VertexType* v0=ep.f->V0(ep.z);
    //            VertexType* v1=ep.f->V1(ep.z);

    //            assert(v0!=v1);

    //            CoordType Pos0=v0->P();
    //            CoordType Pos1=v1->P();

    //            EdgeCoordKey CoordK(std::min(Pos0,Pos1),std::max(Pos0,Pos1));
    //            assert(SplitOps->count(CoordK)>0);
    //            nv.P()=(*SplitOps)[CoordK];
    //        }

    //        vcg::TexCoord2<ScalarType> WedgeInterp(vcg::TexCoord2<ScalarType> &t0,
    //                                               vcg::TexCoord2<ScalarType> &t1)
    //        {
    //            return (vcg::TexCoord2<ScalarType>(0,0));
    //        }

    //        SplitLev(std::map<EdgeCoordKey,CoordType> *_SplitOps){SplitOps=_SplitOps;}
    //    };

    //    class EdgePred
    //    {

    //        std::map<EdgeCoordKey,CoordType> *SplitOps;

    //    public:

    //        bool operator()(vcg::face::Pos<FaceType> ep) const
    //        {
    //            VertexType* v0=ep.f->V0(ep.z);
    //            VertexType* v1=ep.f->V1(ep.z);

    //            assert(v0!=v1);

    //            CoordType Pos0=v0->P();
    //            CoordType Pos1=v1->P();

    //            EdgeCoordKey CoordK(std::min(Pos0,Pos1),std::max(Pos0,Pos1));

    //            return (SplitOps->count(CoordK)>0);
    //        }

    //        EdgePred(std::map<EdgeCoordKey,CoordType> *_SplitOps){SplitOps=_SplitOps;}
    //    };


    //    bool SplitSingularEdgeSeq()
    //    {
    //        std::map<EdgeCoordKey,CoordType> ToBeSplitted;
    //        for (size_t i=0;i<FeatureSeq.size();i++)
    //        {
    //            if (FeatureSeq[i].size()>1)continue;
    //            assert(FeatureSeq[i].size()>0);

    //            //save all edges
    //            CoordType Pos0=FeatureSeq[i][0].V()->P();
    //            CoordType Pos1=FeatureSeq[i][0].VFlip()->P();

    //            EdgeCoordKey SplitOp(std::min(Pos0,Pos1),std::max(Pos0,Pos1));
    //            ToBeSplitted[SplitOp]=(Pos0+Pos1)/2;
    //        }
    //        std::cout<<"Performing "<<ToBeSplitted.size()<< " split ops from Single Edge size Features"<<std::endl;

    //        SplitLev splMd(&ToBeSplitted);
    //        EdgePred eP(&ToBeSplitted);
    //        if (ToBeSplitted.size()==0)return false;
    //        //do the final split
    //        bool done=vcg::tri::RefineE<MeshType,SplitLev,EdgePred>(mesh,splMd,eP);
    //        //update everything
    //        return true;
    //    }

    void CreateEdgeMeshFromEdgeSel(MeshType &Feature)
    {
        std::set<CoordType> SelectedV;
        for (size_t i=0;i<mesh.vert.size();i++)
        {
            if (!mesh.vert[i].IsS())continue;
            SelectedV.insert(mesh.vert[i].P());
        }

        //then create an EdgeMesh
        Feature.Clear();
        for (size_t i=0;i<mesh.face.size();i++)
            for (size_t j=0;j<3;j++)
            {
                if (mesh.face[i].V0(j)>mesh.face[i].V1(j))continue;
                //if (face[i].IsF(j))continue;
                if (!mesh.face[i].IsFaceEdgeS(j))continue;

                CoordType Pos0=mesh.face[i].P0(j);
                CoordType Pos1=mesh.face[i].P1(j);
                vcg::tri::Allocator<MeshType>::AddEdge(Feature,Pos0,Pos1);

            }

        //merge the vertices
        vcg::tri::Clean<MeshType>::RemoveDuplicateVertex(Feature);

        //reselect the vertices
        for (size_t i=0;i<Feature.vert.size();i++)
        {
            if (SelectedV.count(Feature.vert[i].P())==0)continue;
            Feature.vert[i].SetS();
        }
    }

    //    void InitFeatureSeqFromEdgeSel(bool CheckFieldChange,bool checkConcChange)
    //    {
    //        FeatureSeq.clear();
    //        FeatureType.clear();

    //        //first round, get sharp  features to detect cycles
    //        vcg::tri::UpdateSelection<MeshType>::VertexClear(mesh);

    //        InitFaceEdgeSelFromFeatureGeo();

    //        //SelectSplitVertOnSharp(true);
    //        std::cout<<"0 Selected "<<NumVertSel()<<" vertices "<<std::endl;
    //        SelectSplitByMultipleComponents();

    //        if (checkConcChange)
    //            SelectSplitByConcavityChange();

    //        std::cout<<"1 Selected "<<NumVertSel()<<" vertices "<<std::endl;


    //        MeshType Feature;
    //        CreateEdgeMeshFromEdgeSel(Feature);
    //        //std::cout<<"There are Sharp edge mesh:"<<Feature.edge.size()<<std::endl;
    //        //size_t MaxComp=Feature.SplitEdgeConnectedComponents(true);
    //        size_t MaxComp=SplitEdgeConnectedComponents(Feature,true);

    //        InitFeatureSeqFromEdgeSel(MaxComp,Feature);
    //        std::cout<<"There are Sharp sequences:"<<FeatureSeq.size()<<std::endl;

    //        if (!CheckFieldChange)
    //        {
    //            ColorByPartition();
    //            //SetConvexFlatAsConcave(FlatDegree);
    //            return;
    //        }

    //        UpdateEdgeSeqByFieldChange();

    //        //first round, get sharp  features to detect cycles
    //        vcg::tri::UpdateSelection<MeshType>::VertexClear(mesh);
    //        InitFaceEdgeSelFromFeatureGeo();
    //        bool HasSplit=SplitSingularEdgeSeq();
    //        if (HasSplit){
    //            SelectSplitByMultipleComponents();
    //            SelectSplitByConcavityChange();
    //            FeatureSeq.clear();
    //            FeatureType.clear();
    //            Feature.Clear();
    //            CreateEdgeMeshFromEdgeSel(Feature);
    //            //std::cout<<"There are Sharp edge mesh:"<<Feature.edge.size()<<std::endl;
    //            //size_t MaxComp=Feature.SplitEdgeConnectedComponents(true);
    //            size_t MaxComp=SplitEdgeConnectedComponents(Feature,true);

    //            InitFeatureSeqFromEdgeSel(MaxComp,Feature);
    //            InitPreferredDir();
    //            ColorByPartition();
    //            //std::cout<<"There are Sharp sequences:"<<FeatureSeq.size()<<std::endl;
    //        }
    //        //        ColorByPartition();

    //        //SetConvexFlatAsConcave(FlatDegree);
    //    }

    //    void InitFeaturesGeoTableFromEdgeSel()
    //    {
    //        SharpFeaturesGeo.clear();
    //        for (size_t i=0;i<mesh.face.size();i++)
    //            for (size_t j=0;j<3;j++)
    //            {
    //                if (!mesh.face[i].IsFaceEdgeS(j))continue;

    //                CoordType Pos0=mesh.face[i].P0(j);
    //                CoordType Pos1=mesh.face[i].P1(j);

    //                //then add the features
    //                std::pair<CoordType,CoordType> Key(std::min(Pos0,Pos1),std::max(Pos0,Pos1));
    //                SharpFeaturesGeo.insert(Key);
    //            }
    //    }

    //    void InitSharpFeaturesFromCol()
    //    {
    //        for (size_t i=0;i<mesh.face.size();i++)
    //        {
    //            for (size_t j=0;j<mesh.face[i].VN();j++)
    //            {
    //                FaceType *fopp=mesh.face[i].FFp(j);
    //                if (fopp->C()==mesh.face[i].C())
    //                    mesh.face[i].ClearFaceEdgeS(j);
    //                else
    //                    mesh.face[i].SetFaceEdgeS(j);
    //            }
    //        }
    //        InitFeaturesGeoTableFromEdgeSel();
    //    }

    void InitFeatureSequence()
    {
        FeatureSeq.clear();
        FeatureType.clear();

        //first round, get sharp  features to detect cycles
        vcg::tri::UpdateSelection<MeshType>::VertexClear(mesh);

        //SelectSplitVertOnSharp(true);
        // std::cout<<"0 Selected "<<NumVertSel()<<" vertices "<<std::endl;
        SelectSplitByMultipleComponents();

        std::cout<<"0 Selected "<<NumVertSel()<<" vertices "<<std::endl;

        SelectSplitByConcavityChange();

        std::cout<<"1 Selected "<<NumVertSel()<<" vertices "<<std::endl;

        MeshType Feature;
        CreateEdgeMeshFromEdgeSel(Feature);
        size_t MaxComp=MarkEdgeConnectedComponents(Feature,true);

        InitFeatureSeqFromEdgeSel(MaxComp,Feature);
        std::cout<<"There are Sharp sequences:"<<FeatureSeq.size()<<std::endl;

        InitPreferredDir();

        //        ColorByPartition();

        //        //then test connectivity
        //        CheckConnectivity();
        //        if (!CheckFieldChange)
        //        {
        //            ColorByPartition();
        //            //SetConvexFlatAsConcave(FlatDegree);
        //            return;
        //        }

        //        UpdateEdgeSeqByFieldChange();

        //        //first round, get sharp  features to detect cycles
        //        vcg::tri::UpdateSelection<MeshType>::VertexClear(mesh);
        //        InitFaceEdgeSelFromFeatureGeo();
        //        bool HasSplit=SplitSingularEdgeSeq();
        //        if (HasSplit){
        //            SelectSplitByMultipleComponents();
        //            SelectSplitByConcavityChange();
        //            FeatureSeq.clear();
        //            FeatureType.clear();
        //            Feature.Clear();
        //            CreateEdgeMeshFromEdgeSel(Feature);
        //            //std::cout<<"There are Sharp edge mesh:"<<Feature.edge.size()<<std::endl;
        //            //size_t MaxComp=Feature.SplitEdgeConnectedComponents(true);
        //            size_t MaxComp=SplitEdgeConnectedComponents(Feature,true);

        //            InitFeatureSeqFromEdgeSel(MaxComp,Feature);
        //            InitPreferredDir();
        //            ColorByPartition();
        //            //std::cout<<"There are Sharp sequences:"<<FeatureSeq.size()<<std::endl;
        //        }
        //        ColorByPartition();

        //SetConvexFlatAsConcave(FlatDegree);
    }

    bool LoadSharpFeatures(std::string &FeaturePath)
    {
        FILE *f=NULL;
        f=fopen(FeaturePath.c_str(),"rt");
        if(f==NULL) return false;
        int Num=0;
        fscanf(f,"%d\n",&Num);
        std::cout<<"Num "<<Num<<std::endl;
        PerSharpType.clear();
        for (size_t i=0;i<(size_t)Num;i++)
        {
            int FIndex,EIndex;
            int FType;
            fscanf(f,"%d,%d,%d\n",&FType,&FIndex,&EIndex);
            assert(FIndex>=0);
            assert(FIndex<(int)mesh.face.size());
            assert(EIndex>=0);
            assert(EIndex<4);
            CoordType P0=mesh.face[FIndex].P0(EIndex);
            CoordType P1=mesh.face[FIndex].P1(EIndex);
            EdgeCoordKey key(std::min(P0,P1),std::max(P0,P1));
            if (FType==0)
                PerSharpType[key]=ConcaveEdge;
            else
                PerSharpType[key]=ConvexEdge;
            //            if (FType==0)
            //                std::cout<<"Concave"<<std::endl;
            //            else
            //                std::cout<<"Convex"<<std::endl;
            mesh.face[FIndex].SetFaceEdgeS(EIndex);
        }

        InitFeatureSequence();

        ColorByPartition();

        //then test connectivity
        CheckConnectivity();

        CheckFieldFeatureAlignment();

        SortSharpFeatures();

        SplitByDirectionChange();


        return true;
    }

    //    void InitSharpFeatures(ScalarType SharpAngleDegree)/*,
    //                           bool clearTrisAllFeatures=true)*/
    //    {
    //        vcg::tri::UpdateFlags<MeshType>::FaceEdgeSelCrease(mesh,vcg::math::ToRad(SharpAngleDegree));
    ////        if (clearTrisAllFeatures)
    ////        {
    ////            for (size_t i=0;i<mesh.face.size();i++)
    ////            {
    ////                int NumSel=0;
    ////                for (size_t j=0;j<3;j++)
    ////                {
    ////                    if (!mesh.face[i].IsFaceEdgeS(j))continue;
    ////                    NumSel++;
    ////                }
    ////                if (NumSel!=3)continue;
    ////                for (size_t j=0;j<3;j++)
    ////                {
    ////                    if (!IsConcaveEdge(mesh.face[i],j))continue;
    ////                    mesh.face[i].ClearFaceEdgeS(j);

    ////                    FaceType *fOpp=mesh.face[i].FFp(j);
    ////                    int IndexOpp=mesh.face[i].FFi(j);

    ////                    fOpp->ClearFaceEdgeS(IndexOpp);
    ////                }
    ////            }
    ////        }


    //        //        //also add the ones betweek pair of sing
    //        //        for (size_t i=0;i<face.size();i++)
    //        //        {
    //        //            int NumSel=0;
    //        //            for (size_t j=0;j<3;j++)
    //        //            {
    //        //                if (face[i].IsFaceEdgeS(j))continue;
    //        //                bool isSing0=vcg::tri::CrossField<CMesh>::IsSingular(*this,*face[i].V0(j));
    //        //                bool isSing1=vcg::tri::CrossField<CMesh>::IsSingular(*this,*face[i].V1(j));
    //        //                if ((!isSing0)||(!isSing1))continue;

    //        //                face[i].SetFaceEdgeS(j);

    //        //                FaceType *fOpp=face[i].FFp(j);
    //        //                int IndexOpp=face[i].FFi(j);

    //        //                fOpp->SetFaceEdgeS(IndexOpp);
    //        //            }
    //        //        }


    //        InitFeaturesGeoTableFromEdgeSel();
    //    }


    void GetSharpFeatures(std::vector<std::vector<vcg::face::Pos<FaceType> > > &CurrFeatures,
                          const SharpFeatureType CurrFeatureType)
    {
        CurrFeatures.clear();
        for (size_t i=0;i<FeatureSeq.size();i++)
        {
            if (CurrFeatureType!=FeatureType[i])continue;
            CurrFeatures.push_back(FeatureSeq[i]);
        }
    }

    void GetSharpFeatures(std::vector<std::vector<vcg::face::Pos<FaceType> > > &_FeatureSeq,
                          std::vector<SharpFeatureType> &_FeatureType)
    {
        _FeatureSeq=FeatureSeq;
        _FeatureType=FeatureType;
    }

#ifndef NO_TRACING_OPENGL
    void GLDrawSharpFeatures()
    {
        glPushAttrib(GL_ALL_ATTRIB_BITS);
        glDisable(GL_LIGHTING);
        glDepthRange(0,0.9999);
        glLineWidth(5);
        //        glBegin(GL_LINES);
        //        for (size_t i=0;i<FeatureSeq.size();i++)
        //        {
        //            if (ProblematicSeq[i])
        //                vcg::glColor(vcg::Color4b(255,0,0,255));
        //            else
        //            if (FeatureType[i]==ConvexEdge)
        //                vcg::glColor(vcg::Color4b(0,255,0,255));
        //            else
        //                vcg::glColor(vcg::Color4b(0,0,255,255));
        //            for (size_t j=0;j<FeatureSeq[i].size();j++)
        //            {
        //                int VIndex0=vcg::tri::Index(*this,FeatureSeq[i][j].V());
        //                int VIndex1=vcg::tri::Index(*this,FeatureSeq[i][j].VFlip());
        //                CoordType Pos0=vert[VIndex0].P();
        //                CoordType Pos1=vert[VIndex1].P();
        //                vcg::glVertex(Pos0);
        //                vcg::glVertex(Pos1);
        //            }
        //        }
        //        glEnd();
        //        glLineWidth(10);
        //        glDepthRange(0,0.9998);
        //        int done=0;
        glBegin(GL_LINES);
        for (size_t i=0;i<FeatureSeq.size();i++)
        {
            if (FeatureType[i]==ConcaveEdge)
                vcg::glColor(vcg::Color4b(0,0,255,255));
            else
                vcg::glColor(vcg::Color4b(255,0,0,255));

            //vcg::glColor(vcg::Color4b(255,0,0,255));
            //            if (ProblematicSeq[i])
            //                vcg::glColor(vcg::Color4b(255,0,0,255));
            //            else

            //            if (FeatureType[i]!=ConcaveEdge)continue;
            //            done++;
            //            if ((done!=2)&&(done!=3))continue;
            //vcg::glColor(vcg::Color4b::Scatter(FeatureSeq.size(),i));
            for (size_t j=0;j<FeatureSeq[i].size();j++)
            {
                //vcg::glColor(vcg::Color4b::ColorRamp(0,FeatureSeq[i].size(),j));
                int VIndex0=vcg::tri::Index(mesh,FeatureSeq[i][j].V());
                int VIndex1=vcg::tri::Index(mesh,FeatureSeq[i][j].VFlip());
                CoordType Pos0=mesh.vert[VIndex0].P();
                CoordType Pos1=mesh.vert[VIndex1].P();
                vcg::glVertex(Pos0);
                vcg::glVertex(Pos1);
            }
        }
        glEnd();
        glPopAttrib();
    }
#endif

    size_t SelectEndSharpVert(SharpFeatureType SType=NoSharp)
    {
        size_t RetNum=0;
        vcg::tri::UpdateQuality<MeshType>::VertexConstant(mesh,0);
        vcg::tri::UpdateSelection<MeshType>::VertexClear(mesh);
        for (size_t i=0;i<FeatureSeq.size();i++)
        {
            if ((SType!=NoSharp)&&(FeatureType[i]!=SType))continue;
            for (size_t j=0;j<FeatureSeq[i].size();j++)
            {
                int VIndex0=vcg::tri::Index(mesh,FeatureSeq[i][j].V());
                int VIndex1=vcg::tri::Index(mesh,FeatureSeq[i][j].VFlip());
                mesh.vert[VIndex0].Q()+=1;
                mesh.vert[VIndex1].Q()+=1;
            }
        }
        for (size_t i=0;i<mesh.vert.size();i++)
        {
            if (mesh.vert[i].Q()!=1)continue;
            mesh.vert[i].SetS();
            RetNum++;
        }
        return RetNum;
    }

    void SelectSharpConcaveVert()
    {
        vcg::tri::UpdateSelection<MeshType>::VertexClear(mesh);
        for (size_t i=0;i<FeatureSeq.size();i++)
        {
            if (FeatureType[i]!=ConcaveEdge)continue;
            SelectEndSharpVert(FeatureSeq[i],false);
        }
    }

    size_t SelectValence4EndPos()
    {
        size_t RetNum=0;
        vcg::tri::UpdateQuality<MeshType>::VertexConstant(mesh,0);
        vcg::tri::UpdateSelection<MeshType>::VertexClear(mesh);
        for (size_t i=0;i<FeatureSeq.size();i++)
        {
            for (size_t j=0;j<FeatureSeq[i].size();j++)
            {
                int VIndex0=vcg::tri::Index(mesh,FeatureSeq[i][j].V());
                int VIndex1=vcg::tri::Index(mesh,FeatureSeq[i][j].VFlip());
                mesh.vert[VIndex0].Q()+=1;
                mesh.vert[VIndex1].Q()+=1;
            }
        }
        for (size_t i=0;i<mesh.vert.size();i++)
        {
            if (mesh.vert[i].Q()!=4)continue;
            mesh.vert[i].SetS();
            RetNum++;
        }
        return RetNum;
    }

    void SelectSharpVert()
    {
        vcg::tri::UpdateQuality<MeshType>::VertexConstant(mesh,-1);
        vcg::tri::UpdateSelection<MeshType>::VertexClear(mesh);
        for (size_t i=0;i<FeatureSeq.size();i++)
        {
            for (size_t j=0;j<FeatureSeq[i].size();j++)
            {
                int VIndex0=vcg::tri::Index(mesh,FeatureSeq[i][j].V());
                int VIndex1=vcg::tri::Index(mesh,FeatureSeq[i][j].VFlip());
                if (!mesh.vert[VIndex0].IsS())
                {
                    if (mesh.vert[VIndex0].Q()==-1)
                        mesh.vert[VIndex0].Q()=i;
                    if (mesh.vert[VIndex0].Q()!=i)
                        mesh.vert[VIndex0].SetS();
                }
                if (!mesh.vert[VIndex1].IsS())
                {
                    if (mesh.vert[VIndex1].Q()==-1)
                        mesh.vert[VIndex1].Q()=i;
                    if (mesh.vert[VIndex1].Q()!=i)
                        mesh.vert[VIndex1].SetS();
                }
            }
        }

    }

    void SelectEndSharpVert(const std::vector<vcg::face::Pos<FaceType> > &TestSeq,
                            bool clearSel=true)
    {
        vcg::tri::UpdateQuality<MeshType>::VertexConstant(mesh,0);
        if (clearSel)
            vcg::tri::UpdateSelection<MeshType>::VertexClear(mesh);

        for (size_t j=0;j<TestSeq.size();j++)
        {
            int VIndex0=vcg::tri::Index(mesh,TestSeq[j].V());
            int VIndex1=vcg::tri::Index(mesh,TestSeq[j].VFlip());
            mesh.vert[VIndex0].Q()+=1;
            mesh.vert[VIndex1].Q()+=1;
        }
        for (size_t i=0;i<mesh.vert.size();i++)
        {
            if (mesh.vert[i].Q()>=2)continue;
            mesh.vert[i].SetS();
        }
    }

    void SortSharpFeature(std::vector<vcg::face::Pos<FaceType> > &Seq,bool IsLoop)
    {
        if (Seq.size()==1)return;
        vcg::face::Pos<FaceType> StartPos=Seq[0];
        vcg::face::Pos<FaceType> EndPos=Seq[0];
        //        if (IsLoop)
        //            std::cout<<"LOOP"<<std::endl;
        //        else
        //            std::cout<<"NO LOOP"<<std::endl;
        if (!IsLoop)
        {
            bool foundStart=false;
            bool foundEnd=false;
            for (size_t i=0;i<Seq.size();i++)
            {
                if (!foundStart)
                {
                    if (Seq[i].VFlip()->IsS())
                    {
                        StartPos=Seq[i];
                        foundStart=true;
                    }
                    if (Seq[i].V()->IsS())
                    {
                        StartPos=Seq[i];
                        StartPos.FlipV();
                        foundStart=true;
                    }
                }else
                {
                    if (Seq[i].V()->IsS())
                    {
                        EndPos=Seq[i];
                        foundEnd=true;
                    }
                    if (Seq[i].VFlip()->IsS())
                    {
                        EndPos=Seq[i];
                        EndPos.FlipV();
                        foundEnd=true;
                    }
                }
                if (foundStart && foundEnd)break;
            }
            assert(foundStart);
            assert(foundEnd);
        }

        vcg::face::Pos<FaceType> currPos=StartPos;
        std::vector<vcg::face::Pos<FaceType> > NewSeq;
        NewSeq.push_back(currPos);
        do
        {
            VertexType *testV=currPos.V();
            bool found=false;
            for (size_t i=0;i<Seq.size();i++)
            {
                //std::cout<<"i "<<i<<" size 0 "<<Seq.size()<<" size 1 "<<NewSeq.size()<<std::endl;

                if (Seq[i]==currPos)continue;
                vcg::face::Pos<FaceType> oppPos=currPos;
                oppPos.FlipV();
                if (Seq[i]==oppPos)continue;

                if (Seq[i].VFlip()==testV)
                {
                    NewSeq.push_back(Seq[i]);
                    currPos=NewSeq.back();
                    found=true;
                    break;
                }
                if (Seq[i].V()==testV)
                {
                    NewSeq.push_back(Seq[i]);
                    NewSeq.back().FlipV();
                    currPos=NewSeq.back();
                    found=true;
                    break;
                }
            }
            assert(found);
        }while (NewSeq.size()<Seq.size());
        Seq=NewSeq;
    }

    void SortSharpFeatures()
    {

        for (size_t i=0;i<FeatureSeq.size();i++)
        {
            SelectEndSharpVert(FeatureSeq[i]);
            //std::cout<<"index "<<i<<std::endl;
            if (LoopFunctions<MeshType>::IsClosedLoop(FeatureSeq[i]))
            {
                //std::cout<<"Loop "<<std::endl;
                SortSharpFeature(FeatureSeq[i],true);
            }
            else
            {
                //std::cout<<"No Loop "<<std::endl;
                SortSharpFeature(FeatureSeq[i],false);
            }
        }
    }

    //    void SplitAdjacentEdgeSharpFromEdgeSel()
    //    {
    //        vcg::tri::UpdateSelection<MeshType>::VertexClear(mesh);
    //        //InitFaceEdgeSelFromFeatureSeq();

    //        //vcg::tri::UpdateFlags<MeshType>::FaceFauxCrease(*this,math::ToRad(SharpAngleDegree));
    //        //vcg::tri::UpdateFlags<MeshType>::FaceEdgeSelCrease(*this,math::ToRad(SharpAngleDegree));

    //        std::set<std::pair<CoordType,CoordType> > EdgePos;

    //        for (size_t i=0;i<mesh.face.size();i++)
    //        {
    //            for (size_t j=0;j<3;j++)
    //            {
    //                //if (face[i].IsF(j))continue;
    //                if (!mesh.face[i].IsFaceEdgeS(j))continue;
    //                //if (!IsConcaveEdge(face[i],j))continue;
    //                int VIndex0=vcg::tri::Index(mesh,mesh.face[i].V0(j));
    //                int VIndex1=vcg::tri::Index(mesh,mesh.face[i].V1(j));
    //                CoordType P0=mesh.vert[VIndex0].P();
    //                CoordType P1=mesh.vert[VIndex1].P();
    //                mesh.vert[VIndex0].SetS();
    //                mesh.vert[VIndex1].SetS();
    //                EdgePos.insert(std::pair<CoordType,CoordType>(std::min(P0,P1),std::max(P0,P1)));
    //            }
    //        }

    //        //then save the edges to be splitted
    //        std::map<EdgeCoordKey,CoordType> ToBeSplitted;
    //        for (size_t i=0;i<mesh.face.size();i++)
    //        {
    //            //find the number of edges
    //            int Num=0;
    //            for (size_t j=0;j<3;j++)
    //            {
    //                int VIndex0=vcg::tri::Index(mesh,mesh.face[i].V0(j));
    //                int VIndex1=vcg::tri::Index(mesh,mesh.face[i].V1(j));
    //                if ((!mesh.vert[VIndex0].IsS())||(!mesh.vert[VIndex1].IsS()))continue;
    //                CoordType P0=mesh.vert[VIndex0].P();
    //                CoordType P1=mesh.vert[VIndex1].P();
    //                std::pair<CoordType,CoordType> Key(std::min(P0,P1),std::max(P0,P1));
    //                if (EdgePos.count(Key)==1){Num++;continue;}

    //                ToBeSplitted[Key]=(P0+P1)/2;
    //            }
    //            assert(Num<=2);//this should be already solved
    //        }
    //        std::cout<<"Performing "<<ToBeSplitted.size()<< " split ops"<<std::endl;

    //        SplitLev splMd(&ToBeSplitted);
    //        EdgePred eP(&ToBeSplitted);

    //        //do the final split
    //        bool done=vcg::tri::RefineE<MeshType,SplitLev,EdgePred>(mesh,splMd,eP);

    //        mesh.UpdateAttributes();
    //        vcg::tri::CrossField<MeshType>::UpdateSingularByCross(mesh);
    //        //InitFeatureSeqFromEdgeSel();
    //    }

    //    void SplitEdjacentTerminalVertices(bool checkConcavityChange=true)
    //    {
    //        InitFeatureSeqFromEdgeSel(false,checkConcavityChange);

    //        SelectSharpVert();
    //        //vcg::tri::UpdateSelection<MeshType>::VertexClear(*this);
    //        std::set<EdgeCoordKey> NewSelected;
    //        std::map<EdgeCoordKey,CoordType> ToBeSplitted;
    //        for (size_t i=0;i<mesh.face.size();i++)
    //        {
    //            for (size_t j=0;j<3;j++)
    //            {
    //                if (!mesh.face[i].V0(j)->IsS())continue;
    //                if (!mesh.face[i].V1(j)->IsS())continue;

    //                CoordType P0=mesh.face[i].P0(j);
    //                CoordType P1=mesh.face[i].P1(j);
    //                std::pair<CoordType,CoordType> Key(std::min(P0,P1),std::max(P0,P1));

    //                CoordType Mid=(P0+P1)/2;
    //                ToBeSplitted[Key]=Mid;

    //                if (!mesh.face[i].IsFaceEdgeS(j))continue;

    //                std::pair<CoordType,CoordType> Key0(std::min(P0, Mid),std::max(P0, Mid));
    //                std::pair<CoordType,CoordType> Key1(std::min(P1, Mid),std::max(P1, Mid));

    //                NewSelected.insert(Key0);
    //                NewSelected.insert(Key1);
    //            }
    //        }

    //        std::cout<<"Performing "<<ToBeSplitted.size()<< " split ops"<<std::endl;
    //        if (ToBeSplitted.size()==0)return;

    //        SplitLev splMd(&ToBeSplitted);
    //        EdgePred eP(&ToBeSplitted);

    //        //do the final split
    //        bool done=vcg::tri::RefineE<MeshType,SplitLev,EdgePred>(mesh,splMd,eP);

    //        for (size_t i=0;i<mesh.face.size();i++)
    //        {
    //            for (size_t j=0;j<3;j++)
    //            {
    //                if (!mesh.face[i].V0(j)->IsS())continue;
    //                if (!mesh.face[i].V1(j)->IsS())continue;

    //                CoordType P0=mesh.face[i].P0(j);
    //                CoordType P1=mesh.face[i].P1(j);
    //                std::pair<CoordType,CoordType> Key(std::min(P0,P1),std::max(P0,P1));
    //                if (NewSelected.count(Key)==0)continue;
    //                mesh.face[i].SetFaceEdgeS(j);
    //            }
    //        }
    //        InitFeaturesGeoTableFromEdgeSel();
    //    }

    //    void SplitEdgeSharpSharingVerticesFromEdgeSel()
    //    {
    //        int VSize0=mesh.vert.size();
    //        vcg::tri::UpdateColor<MeshType>::PerFaceConstant(mesh);
    //        std::set<EdgeCoordKey> NewSelected;

    //        //vcg::tri::UpdateFlags<MeshType>::FaceFauxCrease(*this,math::ToRad(SharpAngleDegree));
    //        //vcg::tri::UpdateFlags<MeshType>::FaceEdgeSelCrease(*this,math::ToRad(SharpAngleDegree));
    //        //InitFaceEdgeSelFromFeatureSeq();

    //        //save all edges
    //        std::map<EdgeCoordKey,CoordType> ToBeSplitted;
    //        std::set<EdgeCoordKey> MeshEdges;

    //        //safety check
    //        vcg::tri::UpdateSelection<MeshType>::VertexClear(mesh);
    //        for (size_t i=0;i<mesh.face.size();i++)
    //        {
    //            int Num=0;
    //            for (size_t j=0;j<3;j++)
    //            {
    //                CoordType P0=mesh.face[i].P0(j);
    //                CoordType P1=mesh.face[i].P1(j);
    //                std::pair<CoordType,CoordType> Key(std::min(P0,P1),std::max(P0,P1));
    //                MeshEdges.insert(Key);
    //                //if (face[i].IsF(j))continue;
    //                if (!mesh.face[i].IsFaceEdgeS(j))continue;
    //                if (!IsConcaveEdge(mesh.face[i],j))continue;
    //                mesh.face[i].V0(j)->SetS();
    //                mesh.face[i].V1(j)->SetS();
    //                Num++;
    //            }
    //            assert(Num<2);
    //        }

    //        vcg::tri::UpdateFlags<MeshType>::VertexClearV(mesh);
    //        for (size_t i=0;i<mesh.face.size();i++)
    //        {
    //            for (size_t j=0;j<3;j++)
    //            {


    //                FaceType *f0=&mesh.face[i];
    //                FaceType *f1=f0->FFp(j);
    //                int IndexOpp1=f0->FFi(j);
    //                if (f0==f1)continue;

    //                //                if (!f0->IsF(0))continue;
    //                //                if (!f0->IsF(1))continue;
    //                //                if (!f0->IsF(2))continue;
    //                //                if (!f1->IsF(0))continue;
    //                //                if (!f1->IsF(1))continue;
    //                //                if (!f1->IsF(2))continue;

    //                VertexType *Ve0=f0->V0(j);
    //                VertexType *Ve1=f0->V1(j);

    //                VertexType *Vf0=f0->V2(j);
    //                VertexType *Vf1=f1->V2(IndexOpp1);

    //                if (Ve0->IsS())continue;
    //                if (Ve1->IsS())continue;

    //                if (!Vf0->IsS())continue;
    //                if (!Vf1->IsS())continue;

    //                assert(Ve0!=Ve1);
    //                assert(Ve0!=Vf0);
    //                assert(Ve0!=Vf1);
    //                assert(Ve1!=Vf0);
    //                assert(Ve1!=Vf1);
    //                assert(Vf0!=Vf1);

    //                CoordType Pe0=Ve0->P();
    //                CoordType Pe1=Ve1->P();
    //                CoordType Pf0=Vf0->P();
    //                CoordType Pf1=Vf1->P();

    //                EdgeCoordKey SplitOp0(std::min(Pe0,Pf0),std::max(Pe0,Pf0));//edge (j+2)%3 of f0
    //                CoordType Mid0=(SplitOp0.first+SplitOp0.second)/2;
    //                ToBeSplitted[SplitOp0]=Mid0;
    //                if (f0->IsFaceEdgeS((j+2)%3))
    //                {
    //                    EdgeCoordKey E0(std::min(Pe0,Mid0),std::max(Pe0,Mid0));
    //                    EdgeCoordKey E1(std::min(Pf0,Mid0),std::max(Pf0,Mid0));
    //                    NewSelected.insert(E0);
    //                    NewSelected.insert(E1);
    //                }

    //                EdgeCoordKey SplitOp1(std::min(Pe1,Pf0),std::max(Pe1,Pf0));//edge (j+1)%3 of f0
    //                CoordType Mid1=(SplitOp1.first+SplitOp1.second)/2;
    //                ToBeSplitted[SplitOp1]=Mid1;
    //                if (f0->IsFaceEdgeS((j+1)%3))
    //                {
    //                    EdgeCoordKey E0(std::min(Pe1,Mid1),std::max(Pe1,Mid1));
    //                    EdgeCoordKey E1(std::min(Pf0,Mid1),std::max(Pf0,Mid1));
    //                    NewSelected.insert(E0);
    //                    NewSelected.insert(E1);
    //                }

    //                EdgeCoordKey SplitOp2(std::min(Pe0,Pf1),std::max(Pe0,Pf1));//edge (IndexOpp1+1)%3 of f1
    //                CoordType Mid2=(SplitOp2.first+SplitOp2.second)/2;
    //                ToBeSplitted[SplitOp2]=Mid2;
    //                if (f1->IsFaceEdgeS((IndexOpp1+1)%3))
    //                {
    //                    EdgeCoordKey E0(std::min(Pe0,Mid2),std::max(Pe0,Mid2));
    //                    EdgeCoordKey E1(std::min(Pf1,Mid2),std::max(Pf1,Mid2));
    //                    NewSelected.insert(E0);
    //                    NewSelected.insert(E1);
    //                }

    //                EdgeCoordKey SplitOp3(std::min(Pe1,Pf1),std::max(Pe1,Pf1));//edge (IndexOpp1+2)%3 of f1
    //                CoordType Mid3=(SplitOp3.first+SplitOp3.second)/2;
    //                ToBeSplitted[SplitOp3]=Mid3;
    //                if (f1->IsFaceEdgeS((IndexOpp1+2)%3))
    //                {
    //                    EdgeCoordKey E0(std::min(Pe1,Mid3),std::max(Pe1,Mid3));
    //                    EdgeCoordKey E1(std::min(Pf1,Mid3),std::max(Pf1,Mid3));
    //                    NewSelected.insert(E0);
    //                    NewSelected.insert(E1);
    //                }

    //            }
    //        }

    //        std::cout<<"Performing "<<ToBeSplitted.size()<< " split ops"<<std::endl;
    //        if (ToBeSplitted.size()==0)return;

    //        SplitLev splMd(&ToBeSplitted);
    //        EdgePred eP(&ToBeSplitted);

    //        //do the final split
    //        bool done=vcg::tri::RefineE<MeshType,SplitLev,EdgePred>(mesh,splMd,eP);

    //        mesh.UpdateAttributes();
    //        //vcg::tri::CrossField<CMesh>::UpdateSingularByCross(*this);


    //        vcg::tri::UpdateSelection<MeshType>::VertexClear(mesh);
    //        int VSize1=mesh.vert.size();
    //        for (size_t i=VSize0;i<VSize1;i++)
    //            mesh.vert[i].SetS();

    //        for (size_t i=0;i<mesh.face.size();i++)
    //        {
    //            for (size_t j=0;j<3;j++)
    //            {
    //                if (!mesh.face[i].V0(j)->IsS())continue;
    //                if (!mesh.face[i].V1(j)->IsS())continue;

    //                CoordType P0=mesh.face[i].P0(j);
    //                CoordType P1=mesh.face[i].P1(j);
    //                std::pair<CoordType,CoordType> Key(std::min(P0,P1),std::max(P0,P1));
    //                if (NewSelected.count(Key)==0)continue;
    //                mesh.face[i].SetFaceEdgeS(j);
    //            }
    //        }
    //        InitFeaturesGeoTableFromEdgeSel();

    //    }

    //    void RefineInternalFacesStepFromEdgeSel()
    //    {
    //        //vcg::tri::UpdateFlags<MeshType>::FaceEdgeSelCrease(*this,math::ToRad(SharpAngleDegree));
    //        //InitFaceEdgeSelFromFeatureSeq();

    //        std::vector<int> to_refine_face;
    //        for (size_t i=0;i<mesh.face.size();i++)
    //        {
    //            //find the number of edges
    //            int Num=0;
    //            for (size_t j=0;j<3;j++)
    //            {
    //                //if (face[i].IsF(j))continue;
    //                if (!mesh.face[i].IsFaceEdgeS(j))continue;
    //                Num++;
    //            }
    //            if (Num==3)
    //                to_refine_face.push_back(i);
    //        }
    //        if (to_refine_face.size()==0)return;

    //        std::cout<<"WARNING Performing "<<to_refine_face.size()<< " face refinement ops"<<std::endl;
    //        for (size_t j=0;j<to_refine_face.size();j++)
    //        {
    //            int IndexF=to_refine_face[j];
    //            CoordType PD1=mesh.face[IndexF].PD1();
    //            CoordType PD2=mesh.face[IndexF].PD2();
    //            CoordType NewPos=(mesh.face[IndexF].P(0)+mesh.face[IndexF].P(1)+mesh.face[IndexF].P(2))/3;
    //            vcg::tri::Allocator<MeshType>::AddVertex(mesh,NewPos);
    //            VertexType *V0=mesh.face[IndexF].V(0);
    //            VertexType *V1=mesh.face[IndexF].V(1);
    //            VertexType *V2=mesh.face[IndexF].V(2);
    //            VertexType *V3=&mesh.vert.back();
    //            mesh.face[IndexF].V(2)=V3;
    //            vcg::tri::Allocator<MeshType>::AddFace(mesh,V1,V2,V3);
    //            mesh.face.back().PD1()=PD1;
    //            mesh.face.back().PD2()=PD2;
    //            vcg::tri::Allocator<MeshType>::AddFace(mesh,V2,V0,V3);
    //            mesh.face.back().PD1()=PD1;
    //            mesh.face.back().PD2()=PD2;
    //        }
    //        mesh.UpdateAttributes();
    //        //re-update the sequences
    //        //InitFeatureSeqFromEdgeSel();
    //    }

    //    void RefineIfNeeddedFromGeoTable(bool checkConcavityChange=true,bool add_corridor=true)
    //    {
    //        InitFaceEdgeSelFromFeatureGeo();
    //        size_t Size0=mesh.face.size();
    //        RefineInternalFacesStepFromEdgeSel();
    //        size_t Size1=mesh.face.size();
    //        std::cout<<"STEP1 - Triangle with 3 Sharp: Added faces "<<Size1-Size0<<std::endl;

    //        InitFaceEdgeSelFromFeatureGeo();
    //        SplitAdjacentEdgeSharpFromEdgeSel();
    //        size_t Size2=mesh.face.size();
    //        std::cout<<"STEP2 - Triangle with 2 Sharp: Added faces "<<Size2-Size1<<std::endl;

    //        InitFaceEdgeSelFromFeatureGeo();
    //        SplitEdjacentTerminalVertices(checkConcavityChange);
    //        size_t Size3=mesh.face.size();
    //        std::cout<<"STEP 3- Adjacent terminal nodes: Added faces "<<Size3-Size2<<std::endl;

    //        if (add_corridor)
    //        {
    //            InitFaceEdgeSelFromFeatureGeo();
    //            SplitEdgeSharpSharingVerticesFromEdgeSel();
    //            size_t Size4=mesh.face.size();
    //            std::cout<<"STEP 4- Creating corridors "<<Size4-Size3<<std::endl;
    //        }


    //        //        InitFaceEdgeSelFromFeatureGeo();
    //    }

    void CheckAdjacentEdgeSharp()
    {
        //the select the edges
        std::set<std::pair<CoordType,CoordType> > EdgePos;
        for (size_t i=0;i<FeatureSeq.size();i++)
            for (size_t j=0;j<FeatureSeq[i].size();j++)
            {
                int VIndex0=vcg::tri::Index(mesh,FeatureSeq[i][j].V());
                int VIndex1=vcg::tri::Index(mesh,FeatureSeq[i][j].VFlip());
                CoordType P0=mesh.vert[VIndex0].P();
                CoordType P1=mesh.vert[VIndex1].P();
                EdgePos.insert(std::pair<CoordType,CoordType>(std::min(P0,P1),std::max(P0,P1)));
            }
        for (size_t i=0;i<mesh.face.size();i++)
        {
            //find the number of edges
            int Num=0;
            for (size_t j=0;j<3;j++)
            {
                int VIndex0=vcg::tri::Index(mesh,mesh.face[i].V0(j));
                int VIndex1=vcg::tri::Index(mesh,mesh.face[i].V1(j));
                CoordType P0=mesh.vert[VIndex0].P();
                CoordType P1=mesh.vert[VIndex1].P();
                std::pair<CoordType,CoordType> Key(std::min(P0,P1),std::max(P0,P1));
                if (EdgePos.count(Key)==0)continue;
                Num++;
            }
            assert(Num<2);//this should be already solved
        }
    }

    //CHECKING CONNECTIVITY AND FIELD ON SHARP FEATURES
    void CheckSharpEdgesVSVertices()
    {
        vcg::tri::UpdateFlags<MeshType>::VertexClearS(mesh);

        //first check there are no adjacent sharp edhes
        //the select the edges
        std::set<std::pair<CoordType,CoordType> > EdgePos;
        for (size_t i=0;i<FeatureSeq.size();i++)
            for (size_t j=0;j<FeatureSeq[i].size();j++)
            {
                int VIndex0=vcg::tri::Index(mesh,FeatureSeq[i][j].V());
                int VIndex1=vcg::tri::Index(mesh,FeatureSeq[i][j].VFlip());
                mesh.vert[VIndex0].SetS();
                mesh.vert[VIndex1].SetS();
                CoordType P0=mesh.vert[VIndex0].P();
                CoordType P1=mesh.vert[VIndex1].P();
                EdgePos.insert(std::pair<CoordType,CoordType>(std::min(P0,P1),std::max(P0,P1)));
            }

        for (size_t i=0;i<FeatureSeq.size();i++)
            for (size_t j=0;j<FeatureSeq[i].size();j++)
            {
                int VIndex0=vcg::tri::Index(mesh,FeatureSeq[i][j].V());
                int VIndex1=vcg::tri::Index(mesh,FeatureSeq[i][j].VFlip());

                if (mesh.vert[VIndex0].IsS() && mesh.vert[VIndex1].IsS())
                {
                    CoordType P0=mesh.vert[VIndex0].P();
                    CoordType P1=mesh.vert[VIndex1].P();
                    std::pair<CoordType,CoordType> Key(std::min(P0,P1),std::max(P0,P1));
                    assert(EdgePos.count(Key)==1);
                }
            }
    }

    //bool colored;

    //    void CheckCorridorBetweenSharp()
    //    {
    //        colored=false;

    //        //safety check
    //        vcg::tri::UpdateSelection<MeshType>::VertexClear(*this);
    //        for (size_t i=0;i<FeatureSeq.size();i++)
    //            for (size_t j=0;j<FeatureSeq[i].size();j++)
    //            {
    //                int VIndex0=vcg::tri::Index(*this,FeatureSeq[i][j].V());
    //                int VIndex1=vcg::tri::Index(*this,FeatureSeq[i][j].VFlip());
    //                if (!IsConcaveEdge(*FeatureSeq[i][j].F(),FeatureSeq[i][j].E()))continue;
    //                mesh.vert[VIndex0].SetS();
    //                mesh.vert[VIndex1].SetS();
    //            }


    //        vcg::tri::UpdateFlags<MeshType>::VertexClearV(*this);
    //        for (size_t i=0;i<mesh.face.size();i++)
    //        {
    //            for (size_t j=0;j<3;j++)
    //            {
    //                FaceType *f0=&mesh.face[i];
    //                FaceType *f1=f0->FFp(j);
    //                int IndexOpp1=f0->FFi(j);
    //                if (f0==f1)continue;

    //                VertexType *Ve0=f0->V0(j);
    //                VertexType *Ve1=f0->V1(j);

    //                VertexType *Vf0=f0->V2(j);
    //                VertexType *Vf1=f1->V2(IndexOpp1);

    //                if (Ve0->IsS())continue;
    //                if (Ve1->IsS())continue;

    //                if (!Vf0->IsS())continue;
    //                if (!Vf1->IsS())continue;

    //                if (colored)continue;
    //                colored=true;
    //                f0->C()=vcg::Color4b::Red;
    //                f1->C()=vcg::Color4b::Red;
    //                std::cout<<"WARNING NON CORRIDOR BETWEEN FEATURES"<<std::endl;
    //                assert(0);

    //            }
    //        }
    //    }

    void CheckNoFaceAllSharp()
    {
        vcg::tri::UpdateFlags<MeshType>::FaceClearFaceEdgeS(mesh);
        for (size_t i=0;i<FeatureSeq.size();i++)
            for (size_t j=0;j<FeatureSeq[i].size();j++)
            {
                vcg::face::Pos<FaceType> CurrPos=FeatureSeq[i][j];
                CurrPos.F()->SetFaceEdgeS(CurrPos.E());
                CurrPos.FlipF();
                CurrPos.F()->SetFaceEdgeS(CurrPos.E());
            }
        for (size_t i=0;i<mesh.face.size();i++)
        {
            int NumSel=0;
            for (size_t j=0;j<3;j++)
                if (mesh.face[i].IsFaceEdgeS(j))NumSel++;
            assert(NumSel<3);
        }
    }

    void CheckConnectivity()
    {
        CheckSharpEdgesVSVertices();
        //CheckCorridorBetweenSharp();
        CheckNoFaceAllSharp();
    }

    //    void RemoveSmallFeatures(size_t MinSize=2)
    //    {
    //        SelectEndSharpVert();
    //        std::vector<bool> HasEndPoint(FeatureSeq.size(),false);
    //        for (size_t i=0;i<FeatureSeq.size();i++)
    //            for (size_t j=0;j<FeatureSeq[i].size();j++)
    //            {
    //                VertexType *v0=FeatureSeq[i][j].V();
    //                VertexType *v1=FeatureSeq[i][j].VFlip();
    //                if (v0->IsS()){HasEndPoint[i]=true;break;}
    //                if (v1->IsS()){HasEndPoint[i]=true;break;}
    //            }

    //        std::vector<std::vector<vcg::face::Pos<FaceType> > > FeatureSeqTemp;
    //        std::vector<SharpFeatureType> FeatureTypeTemp;
    //        std::vector<bool> ProblematicSeqTemp;

    //        size_t Num0=FeatureSeq.size();
    //        for (size_t i=0;i<FeatureSeq.size();i++)
    //        {
    //            if ((HasEndPoint[i])&&(FeatureSeq[i].size()<=MinSize))continue;

    //            FeatureSeqTemp.push_back(FeatureSeq[i]);
    //            FeatureTypeTemp.push_back(FeatureType[i]);
    //            ProblematicSeqTemp.push_back(ProblematicSeq[i]);
    //        }

    //        FeatureSeq=FeatureSeqTemp;
    //        FeatureType=FeatureTypeTemp;
    //        ProblematicSeq=ProblematicSeqTemp;

    //        size_t Num1=FeatureSeq.size();

    //        vcg::tri::UpdateFlags<MeshType>::FaceClearFaceEdgeS(mesh);
    //        for (size_t i=0;i<FeatureSeq.size();i++)
    //            for (size_t j=0;j<FeatureSeq[i].size();j++)
    //            {
    //                vcg::face::Pos<FaceType> currPos=FeatureSeq[i][j];
    //                FaceType *f0=currPos.F();
    //                int E0=currPos.E();
    //                currPos.FlipF();
    //                FaceType *f1=currPos.F();
    //                int E1=currPos.E();
    //                f0->SetFaceEdgeS(E0);
    //                f1->SetFaceEdgeS(E1);
    //            }
    //        InitFeaturesGeoTableFromEdgeSel();
    //        std::cout<<std::endl<<std::endl<<"*** Removed "<<Num0-Num1<<" sharp sequences ***"<<std::endl;
    //    }

    bool CorrectAlignment(FaceType *F,const int &IndexE,ScalarType Tolerance=0.01)
    {
        CoordType Pos0=F->P0(IndexE);
        CoordType Pos1=F->P1(IndexE);
        CoordType Dir=Pos1-Pos0;
        Dir.Normalize();
        CoordType PD1=F->PD1();
        PD1.Normalize();
        if (fabs(PD1*Dir)<Tolerance)return false;
        if (fabs(PD1*Dir)>(1-Tolerance))return false;
        F->C()=vcg::Color4b(255,0,0,255);
        //F->PD1()=Dir;
        //F->PD2()=Dir^F->N();
        return true;
    }

    void CheckFieldFeatureAlignment()
    {
        int Modified=0;
        for (size_t i=0;i<FeatureSeq.size();i++)
            for (size_t j=0;j<FeatureSeq[i].size();j++)
            {
                vcg::face::Pos<FaceType> CurrPos=FeatureSeq[i][j];
                bool has_modified=CorrectAlignment(CurrPos.F(),CurrPos.E());
                if (has_modified) Modified++;
                CurrPos.FlipF();
                has_modified=CorrectAlignment(CurrPos.F(),CurrPos.E());
                if (has_modified) Modified++;
            }
        if (Modified==0)return;
        std::cout<<"WARNING modified field "<<Modified<<std::endl;
    }

    ScalarType AvgAngle(size_t IndexSeq)
    {
        ScalarType SumAngle=0;
        for (size_t j=0;j<FeatureSeq[IndexSeq].size();j++)
        {
            ScalarType angle = vcg::face::DihedralAngleRad(*FeatureSeq[IndexSeq][j].F(),FeatureSeq[IndexSeq][j].E());
            angle=vcg::math::ToDeg(angle);
            //std::cout<<"Angle test"<<angle<<std::endl;
            SumAngle+=angle;
        }
        return (SumAngle/FeatureSeq[IndexSeq].size());
    }

    //    void SetConvexFlatAsConcave(ScalarType minAngle)
    //    {
    //        // get average angle
    //        for (size_t i=0;i<FeatureSeq.size();i++)
    //        {
    //            if (FeatureType[i]!=ConvexEdge)continue;
    //            if (AvgAngle(i)>minAngle)continue;
    //            //std::cout<<"Angle "<<AvgAngle(i)<<std::endl;
    //            FeatureType[i]=ConcaveEdge;
    //        }
    //    }
};
#endif
