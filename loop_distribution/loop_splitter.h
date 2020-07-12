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

#ifndef LOOP_SPLITTER_H
#define LOOP_SPLITTER_H

#include "tracing_field/sharp_feature_manager.h"
#include "tracing_field/remesh/edge_mesh_type.h"

template <class MeshType>
class LoopSplitter
{


private:

    typedef typename MeshType::CoordType    CoordType;
    typedef typename MeshType::VertexType   VertexType;
    typedef typename MeshType::FaceType     FaceType;
    typedef typename MeshType::ScalarType   ScalarType;
    typedef Path<ScalarType> PathType;

    AnisotropicGraph<MeshType> &anigraph;
    SharpFeaturesManager<MeshType> &sharp;

    typedef std::pair<size_t,size_t> FaceEdgePair;
    typedef std::pair<CoordType,CoordType> CoordPair;
    MyEMesh EdgeM;
    std::vector<std::vector<FaceEdgePair> > SideConcaveFeatures;
    std::vector<std::vector<FaceEdgePair> > SideConvexFeatures;

public:

    enum LoopType {Concave,ConcaveNonSampled,ConvexFeature,Regular,Undetermined};

    struct LoopSeq{
        std::vector<FaceEdgePair> IndexFE;
        std::vector<bool> Sharp;
        bool IsLoop;
        LoopType LType;
        bool CrossProperly;
        //        inline bool operator <(const LoopSeq &LSeq)const
        //        {
        //            assert(LType!=Undetermined);
        //            assert(LSeq.LType!=Undetermined);
        //            if (LType!=LSeq.LType)
        //            {
        //                if (LType==Concave)return true;
        //                if (LSeq.LType==Concave)return false;
        //                if (LType==ConcaveNonSampled)return true;
        //                if (LSeq.LType==ConcaveNonSampled)return false;
        //                if (LType==ConvexFeature)return true;
        //                if (LSeq.LType==ConvexFeature)return false;
        //                assert(0);
        //            }
        //            if (IsLoop!=LSeq.IsLoop)
        //            {
        //                if (IsLoop)return true;
        //                if (LSeq.IsLoop)return false;
        //                assert(0);
        //            }
        //        }

        LoopSeq()
        {
            IsLoop=false;
            LType=Undetermined;
            CrossProperly=true;
        }
    };

    std::vector<LoopSeq> LoopSequences;

private:

    void GetEdgeFaceOnFeatures(const std::vector<size_t> &IndexF,
                               std::vector<FaceEdgePair> &OnFeaturePos)
    {
        OnFeaturePos.clear();
        //std::cout<<"Size "<<IndexF.size()<<std::endl;
        for (size_t j=0;j<IndexF.size();j++)
        {
            int sharpE=-1;
            size_t currF=IndexF[j];
            for (size_t e=0;e<3;e++)
            {
                if (!anigraph.Mesh().face[currF].IsFaceEdgeS(e))continue;
                assert(sharpE==-1);//Only one feature per face
                sharpE=e;
            }
            if (sharpE==-1)continue;//face which not insist on sharp feature
            OnFeaturePos.push_back(FaceEdgePair(currF,sharpE));
        }
        assert(OnFeaturePos.size()>0);//there should be at least one
    }

    void InitSideSequences()
    {
        SideConcaveFeatures.clear();
        SideConvexFeatures.clear();

        assert(sharp.IndexF0.size()==sharp.IndexF1.size());
        assert(sharp.FeatureType.size()==sharp.IndexF0.size());

        sharp.InitFaceEdgeSelFromEdgeSeq();
        for (size_t i=0;i<sharp.IndexF0.size();i++)
        {
            if (sharp.FeatureType[i]==ConcaveEdge)
            {
                SideConcaveFeatures.resize(SideConcaveFeatures.size()+1);
                GetEdgeFaceOnFeatures(sharp.IndexF0[i],SideConcaveFeatures.back());

                SideConcaveFeatures.resize(SideConcaveFeatures.size()+1);
                GetEdgeFaceOnFeatures(sharp.IndexF1[i],SideConcaveFeatures.back());
            }
            if (sharp.FeatureType[i]==ConvexEdge)
            {
                SideConvexFeatures.resize(SideConvexFeatures.size()+1);
                GetEdgeFaceOnFeatures(sharp.IndexF0[i],SideConvexFeatures.back());

                SideConvexFeatures.resize(SideConvexFeatures.size()+1);
                GetEdgeFaceOnFeatures(sharp.IndexF1[i],SideConvexFeatures.back());
            }
        }
    }


    struct FaceEdge
    {
        CoordType P0;
        CoordType P1;
        size_t IndexF;

        inline bool operator ==(const FaceEdge &Fedge)const
        {
            return ((P0==Fedge.P0)&&(P1==Fedge.P1)&&(IndexF==Fedge.IndexF));
        }

        inline bool operator <(const FaceEdge &Fedge)const
        {
            if (P0!=Fedge.P0)return (P0<Fedge.P0);
            if (P1!=Fedge.P1)return (P1<Fedge.P1);
            return (IndexF<Fedge.IndexF);
        }

        FaceEdge(const CoordType &_P0,
                 const CoordType &_P1,
                 const size_t &_IndexF)
        {
            P0=std::min(_P0,_P1);
            P1=std::max(_P0,_P1);
            IndexF=_IndexF;
        }

    };

    void ColorByPathIndex(const std::vector<PathType> &ChoosenPaths)
    {
        vcg::tri::UpdateColor<MeshType>::PerFaceConstant(anigraph.Mesh(),vcg::Color4b::White);
        for (size_t i=0;i<anigraph.Mesh().face.size();i++)
        {
            if (anigraph.Mesh().face[i].IndexFeature==-1)continue;
            int IndexPath=anigraph.Mesh().face[i].FeaturePath;
            vcg::Color4b col;
            if (!anigraph.Mesh().face[i].Concave)
            {
                assert(IndexPath==-1);
                col=vcg::Color4b::Green;
            }
            else{
                //in case not sampled
                if (IndexPath<0)
                {
                    col=vcg::Color4b::Red;
                }else
                {
                    col=vcg::Color4b::Scatter(ChoosenPaths.size(),IndexPath);
                }
            }
            anigraph.Mesh().face[i].C()=col;
        }
    }

    void UpdatePathForSequenceAfterSplit()
    {

        for (size_t i=0;i<SideConcaveFeatures.size();i++)
            SideConcaveFeatures[i].clear();
        for (size_t i=0;i<SideConvexFeatures.size();i++)
            SideConvexFeatures[i].clear();

        //update the feature
        for (size_t i=0;i<anigraph.Mesh().face.size();i++)
        {
            int currFeature=anigraph.Mesh().face[i].IndexFeature;
            if (currFeature==-1)continue;

            //find the next adjacent (check what it should be next feature)
            int nextFeature=currFeature+1;
            if (currFeature%2==1)
                nextFeature=currFeature-1;

            bool IsConcave=anigraph.Mesh().face[i].Concave;
            int foundEdge=-1;
            for (size_t j=0;j<3;j++)
            {
                FaceType *nextF=anigraph.Mesh().face[i].FFp(j);
                int testFeature=nextF->IndexFeature;
                bool IsConcaveTest=nextF->Concave;
                if ((testFeature==nextFeature)&&(IsConcaveTest==IsConcave))
                {
                    assert(foundEdge==-1);
                    foundEdge=j;
                }
            }
            if (foundEdge==-1)continue;

            anigraph.Mesh().face[i].IndexE=foundEdge;
            if (anigraph.Mesh().face[i].Concave)
            {
                assert(currFeature>=0);
                assert(currFeature<(int)SideConcaveFeatures.size());
                assert(nextFeature>=0);
                assert(nextFeature<(int)SideConcaveFeatures.size());
                SideConcaveFeatures[currFeature].push_back(FaceEdgePair(i,foundEdge));
            }
            else
            {
                assert(currFeature>=0);
                assert(currFeature<(int)SideConvexFeatures.size());
                assert(nextFeature>=0);
                assert(nextFeature<(int)SideConvexFeatures.size());
                SideConvexFeatures[currFeature].push_back(FaceEdgePair(i,foundEdge));
            }
        }

        for (size_t i=0;i<anigraph.Mesh().face.size();i++)
        {
            anigraph.Mesh().face[i].Concave=false;
            anigraph.Mesh().face[i].IndexFeature=-1;
        }

        for (size_t i=0;i<SideConcaveFeatures.size();i++)
            for (size_t j=0;j<SideConcaveFeatures[i].size();j++)
            {
                size_t IndexF=SideConcaveFeatures[i][j].first;
                anigraph.Mesh().face[IndexF].IndexFeature=i;
                anigraph.Mesh().face[IndexF].Concave=true;
            }

        for (size_t i=0;i<SideConvexFeatures.size();i++)
            for (size_t j=0;j<SideConvexFeatures[i].size();j++)
            {
                size_t IndexF=SideConvexFeatures[i][j].first;
                anigraph.Mesh().face[IndexF].IndexFeature=i;
                anigraph.Mesh().face[IndexF].Concave=false;
            }
        for (size_t i=0;i<anigraph.Mesh().face.size();i++)
        {
            if (anigraph.Mesh().face[i].IndexFeature==-1)
                anigraph.Mesh().face[i].FeaturePath=-1;
        }
    }

    void InitFacePathForSequence(const std::vector<PathType> &ChoosenPaths,
                                 const std::vector<bool> &IsLoop)
    {
        //        FacePath.clear();
        //        FacePath.resize(anigraph.Mesh().face.size(),-1);

        for (size_t i=0;i<anigraph.Mesh().face.size();i++)
        {
            anigraph.Mesh().face[i].Concave=false;
            anigraph.Mesh().face[i].FeaturePath=-1;
            anigraph.Mesh().face[i].IndexFeature=-1;
        }

        for (size_t i=0;i<SideConcaveFeatures.size();i++)
        {
            for (size_t j=0;j<SideConcaveFeatures[i].size();j++)
            {
                size_t IndexF=SideConcaveFeatures[i][j].first;
                anigraph.Mesh().face[IndexF].Concave=true;
                assert(anigraph.Mesh().face[IndexF].IndexFeature==-1);
                anigraph.Mesh().face[IndexF].IndexFeature=i;
            }
        }

        for (size_t i=0;i<SideConvexFeatures.size();i++)
        {
            for (size_t j=0;j<SideConvexFeatures[i].size();j++)
            {
                size_t IndexF=SideConvexFeatures[i][j].first;
                anigraph.Mesh().face[IndexF].Concave=false;
                assert(anigraph.Mesh().face[IndexF].IndexFeature==-1);
                anigraph.Mesh().face[IndexF].IndexFeature=i;
            }
        }

        //first initialize a map for loop edge/face
        std::map<FaceEdge,size_t> Fedges;

        for (size_t i=0;i<ChoosenPaths.size();i++)
        {
            std::vector<CoordType> PathPos;
            anigraph.PathPos(ChoosenPaths[i],PathPos);
            std::vector<std::pair<size_t,size_t> > FaceDirTemp;
            anigraph.GetFacesDir(ChoosenPaths[i],FaceDirTemp,IsLoop[i]);

            //size_t limit=FaceDirTemp.size();
            if (IsLoop[i])
            {
                assert(PathPos.size()==(FaceDirTemp.size()));
            }
            else
            {
                assert(PathPos.size()==(FaceDirTemp.size()+1));
                //limit=FaceDirTemp.size()-1;
            }

            for (size_t j=0;j<FaceDirTemp.size();j++)
            {
                size_t IndexF=FaceDirTemp[j].first;
                assert(IndexF>=0);
                assert(IndexF<anigraph.Mesh().face.size());
                CoordType P0=PathPos[j];
                CoordType P1=PathPos[(j+1)%PathPos.size()];
                if (P0==P1)continue;//snapped on the same position
                FaceEdge key(P0,P1,IndexF);
                //                if (Fedges.count(key)>0)
                //                {
                //                    std::cout<<"Pos0 "<<P0.X()<<","<<P0.Y()<<","<<P0.Z()<<std::endl;
                //                    std::cout<<"Pos1 "<<P1.X()<<","<<P1.Y()<<","<<P1.Z()<<std::endl;
                //                    std::cout<<"Path0 "<<Fedges[key]<<std::endl;
                //                    std::cout<<"Path1 "<<i<<std::endl;
                //                    assert(0);
                //                }
                assert(Fedges.count(key)==0);
                Fedges[key]=i;
            }
        }

        std::vector<int> PathConcaveIndex;
        PathConcaveIndex.resize(SideConcaveFeatures.size(),-2);

        //for each sequence
        for (size_t i=0;i<SideConcaveFeatures.size();i++)
        {
            //for each edge
            for (size_t j=0;j<SideConcaveFeatures[i].size();j++)
            {
                size_t IndexF=SideConcaveFeatures[i][j].first;
                size_t IndexE=SideConcaveFeatures[i][j].second;
                assert(IndexF>=0);
                assert(IndexF<anigraph.Mesh().face.size());
                assert(IndexE>=0);
                assert(IndexE<3);
                CoordType P0=anigraph.Mesh().face[IndexF].P0(IndexE);
                CoordType P1=anigraph.Mesh().face[IndexF].P1(IndexE);
                FaceEdge key(P0,P1,IndexF);
                if (Fedges.count(key)==0)
                    PathConcaveIndex[i]=-1;
                else
                {
                    assert(PathConcaveIndex[i]!=-1);
                    size_t currP=Fedges[key];
                    if (PathConcaveIndex[i]==-2)
                        PathConcaveIndex[i]=currP;
                    else{
                        assert(PathConcaveIndex[i]==(int)currP);
                    }
                }
            }
        }

        //then set the path on the quality the first N are the regular
        //features, then add the others
        for (size_t i=0;i<SideConcaveFeatures.size();i++)
        {
            assert(PathConcaveIndex[i]!=-2);
            //for each edge
            for (size_t j=0;j<SideConcaveFeatures[i].size();j++)
            {
                size_t IndexF=SideConcaveFeatures[i][j].first;
                //if (PathConcaveIndex[i]>=0)
                anigraph.Mesh().face[IndexF].FeaturePath=PathConcaveIndex[i];
                //anigraph.Mesh().face[IndexF].Q()=PathConcaveIndex[i];
            }
        }
        //then color the faces
        //ColorByPathQuality(ChoosenPaths);
    }


    void AllocateAdditionalParameters(MeshType &m)
    {
        if (!vcg::tri::HasPerVertexAttribute(m,std::string("OriginalV")))
            vcg::tri::Allocator<MeshType>:: template AddPerVertexAttribute<bool> (m,std::string("OriginalV"));
    }

    struct SharpLoopInfo
    {
        int IndexF0;
        int IndexF1;
        int IndexL0;
        int IndexL1;

        SharpLoopInfo(const int &_IndexF0,
                      const int &_IndexF1)
        {
            IndexF0=_IndexF0;
            IndexF1=_IndexF1;
            assert(IndexF1=IndexF0+1);
            IndexL0=-1;
            IndexL1=-1;
        }
    };

    void SplitMesh(const std::vector<PathType> &ChoosenPaths,
                   const std::vector<bool> &IsLoop)
    {
        EdgeM.Clear();

        AllocateAdditionalParameters(anigraph.Mesh());

        std::cout<<"getting Edge Mesh"<<std::endl;
        anigraph.GetEdgeMesh(ChoosenPaths,IsLoop,EdgeM,anigraph.Mesh());

        //remove the sharp features
        EdgeM.RemoveEdgeFromMesh(anigraph.Mesh());

        //        vcg::tri::io::ExporterPLY<MyEMesh>::Save(EdgeM,"./testE1.ply",vcg::tri::io::Mask::IOM_EDGEINDEX);
        EdgeSplitter< MeshType,MyEMesh > esplit(anigraph.Mesh(),EdgeM);

        std::cout<<"Splitting Internal"<<std::endl;
        esplit.SplitInternal();

        std::cout<<"Splitting The rest"<<std::endl;
        //EdgeSplitter<MeshType,MyEMesh> EdgeSplit(anigraph.Mesh(),EdgeM);
        //EdgeSplit.Split(false);

        esplit.Split(false);
        std::cout<<"Splitting Step End"<<std::endl;

        //        vcg::tri::io::ExporterPLY<MyEMesh>::Save(EdgeM,"./testE.ply",vcg::tri::io::Mask::IOM_EDGEINDEX);

    }

    void RetrievePath(VertexType *v0,
                      VertexType *v1,
                      std::vector<CoordType> &PosPath)
    {
        PosPath.clear();
        std::vector<VertexType*> Seed(1,v0);
        EuclideanDistance<MeshType> dd;
        typename MeshType::template PerVertexAttributeHandle<VertexType*> parentV =
                vcg::tri::Allocator<MeshType>::template GetPerVertexAttribute<VertexType*> (anigraph.Mesh());
        ScalarType MaxD=anigraph.Mesh().bbox.Diag();
        vcg::tri::Geodesic<MeshType>::PerVertexDijkstraCompute(anigraph.Mesh(),Seed,dd,MaxD,NULL,NULL,&parentV,false,v1);
        PosPath.push_back(v1->P());
        VertexType *CurrV=v1;
        do
        {
            CurrV=parentV[CurrV];
            PosPath.push_back(CurrV->P());
        }while (CurrV!=v0);
        std::reverse(PosPath.begin(),PosPath.end());
    }

    void InitLoopSequences(const std::vector<bool> &CrossOK)
    {
        //vcg::tri::io::ExporterPLY<MyEMesh>::Save(EdgeM,"./testE1.ply",vcg::tri::io::Mask::IOM_EDGEINDEX);
        //vcg::tri::io::ExporterPLY<MeshType>::Save(anigraph.Mesh(),"./testM1.ply");


        LoopSequences.clear();

        //create a map with mesh edge
        std::map<CoordPair,FaceEdgePair> EdgeMap;
        for (size_t i=0;i<anigraph.Mesh().face.size();i++)
            for (size_t j=0;j<3;j++)
            {
                CoordType Pos0=anigraph.Mesh().face[i].P0(j);
                CoordType Pos1=anigraph.Mesh().face[i].P1(j);
                CoordPair key(std::min(Pos0,Pos1),std::max(Pos0,Pos1));
                EdgeMap[key]=FaceEdgePair(i,j);
            }
        vcg::GridStaticPtr<VertexType,ScalarType> Gr;
        Gr.Set(anigraph.Mesh().vert.begin(),anigraph.Mesh().vert.end());

        //set one loop after the other
        int MaxPath=0;
        for (size_t i=0;i<EdgeM.edge.size();i++)
            MaxPath=std::max(MaxPath,EdgeM.edge[i].PathIndex);

        LoopSequences.resize(MaxPath+1);
        for (size_t i=0;i<EdgeM.edge.size();i++)
        {
            size_t currI=EdgeM.edge[i].PathIndex;
            assert(currI<LoopSequences.size());
            assert(currI>=0);
            //then  retrieve the face
            CoordType Pos0=EdgeM.edge[i].P(0);
            ScalarType minD;
            VertexType *v0=vcg::tri::GetClosestVertex(anigraph.Mesh(),Gr,Pos0,anigraph.Mesh().bbox.Diag(),minD);
            CoordType Pos0t=v0->P();
            CoordType Pos1=EdgeM.edge[i].P(1);
            VertexType *v1=vcg::tri::GetClosestVertex(anigraph.Mesh(),Gr,Pos1,anigraph.Mesh().bbox.Diag(),minD);
            CoordType Pos1t=v1->P();
            //            CoordPair key(std::min(Pos0,Pos1),std::max(Pos0,Pos1));
            //            std::cout<<"C "<<EdgeMap.count(key)<<std::endl;
            //            std::cout<<"Pos0 "<<Pos0.X()<<","<<Pos0.Y()<<","<<Pos0.Z()<<","<<std::endl;
            //            std::cout<<"Pos0t "<<Pos0t.X()<<","<<Pos0t.Y()<<","<<Pos0t.Z()<<","<<std::endl;
            //            std::cout<<"Pos1 "<<Pos1.X()<<","<<Pos1.Y()<<","<<Pos1.Z()<<","<<std::endl;
            //            std::cout<<"Pos1t "<<Pos1t.X()<<","<<Pos1t.Y()<<","<<Pos1t.Z()<<","<<std::endl;
            CoordPair key1(std::min(Pos0t,Pos1t),std::max(Pos0t,Pos1t));
            //std::cout<<"C1 "<<EdgeMap.count(key1)<<std::endl;
            if(EdgeMap.count(key1)==0)
            {
                //still possible in some case
                std::vector<CoordType> PosPath;
                RetrievePath(v0,v1,PosPath);
                assert(PosPath.size()>2);
                for (size_t j=0;j<PosPath.size()-1;j++)
                {
                    CoordType P0=PosPath[j];
                    CoordType P1=PosPath[j+1];
                    CoordPair key2(std::min(P0,P1),std::max(P0,P1));
                    assert(EdgeMap.count(key2)>0);
                    FaceEdgePair currFE=EdgeMap[key2];
                    LoopSequences[currI].IndexFE.push_back(currFE);
                    LoopSequences[currI].Sharp.push_back(false);
                }
                //                MyEMesh TestM;
                //                vcg::tri::Allocator<MyEMesh>::AddEdge(TestM,Pos0t,Pos1t);
                //                vcg::tri::io::ExporterPLY<MeshType>::Save(anigraph.Mesh(),"./SplittedM.ply");
                //                vcg::tri::io::ExporterPLY<MyEMesh>::Save(TestM,"./testSingleE.ply",vcg::tri::io::Mask::IOM_EDGEINDEX);
                //                assert(0);
            }else
            {
                FaceEdgePair currFE=EdgeMap[key1];
                LoopSequences[currI].IndexFE.push_back(currFE);
                LoopSequences[currI].Sharp.push_back(false);
            }
        }

        //set this one as regular
        for (size_t i=0;i<LoopSequences.size();i++)
        {
            LoopSequences[i].LType=Regular;
            LoopSequences[i].CrossProperly=CrossOK[i];
            //at least one edge per sequence
            //assert(LoopSequences[i].IndexFE.size()>0);
        }

        std::vector<bool> IsSampled(SideConcaveFeatures.size(),false);
        for (size_t i=0;i<anigraph.Mesh().face.size();i++)
        {
            int IndexFeature=anigraph.Mesh().face[i].IndexFeature;
            int FeaturePath=anigraph.Mesh().face[i].FeaturePath;
            int IndexE=anigraph.Mesh().face[i].IndexE;
            bool IsConcave=anigraph.Mesh().face[i].Concave;
            if (!IsConcave)continue;
            if (FeaturePath==-1)continue;
            //sampled path
            assert(IndexFeature!=-1);
            assert(IndexFeature>=0);
            assert(IndexFeature<(int)IsSampled.size());
            assert(FeaturePath>=0);
            assert(FeaturePath<(int)LoopSequences.size());
            LoopSequences[FeaturePath].IndexFE.push_back(FaceEdgePair(i,IndexE));
            LoopSequences[FeaturePath].Sharp.push_back(true);
            LoopSequences[FeaturePath].LType=Concave;
            IsSampled[IndexFeature]=true;
        }

        //then the concave ones that have not been sampled
        for (size_t i=0;i<SideConcaveFeatures.size();i++)
        {
            if (IsSampled[i])continue;//already inserted
            LoopSequences.resize(LoopSequences.size()+1);
            LoopSequences.back().LType=ConcaveNonSampled;
            LoopSequences.back().IndexFE=SideConcaveFeatures[i];
            LoopSequences.back().Sharp.resize(LoopSequences.back().IndexFE.size(),true);
        }

        //set the convex ones
        for (size_t i=0;i<SideConvexFeatures.size();i++)
        {
            LoopSequences.resize(LoopSequences.size()+1);
            LoopSequences.back().LType=ConvexFeature;
            for (size_t j=0;j<SideConvexFeatures[i].size();j++)
            {
                int CurrF=SideConvexFeatures[i][j].first;
                int CurrE=SideConvexFeatures[i][j].second;
                LoopSequences.back().IndexFE.push_back(FaceEdgePair(CurrF,CurrE));
                LoopSequences.back().Sharp.push_back(true);
            }
        }

        //then check the loops and resort to put complete concave sharp at the beginning
        std::vector<bool> Added(LoopSequences.size(),false);
        std::vector<LoopSeq> LoopSwap;
        for (size_t i=0;i<LoopSequences.size();i++)
        {
            LoopSequences[i].IsLoop=LoopFunctions<MeshType>::IsClosedLoop(anigraph.Mesh(),LoopSequences[i].IndexFE);
            //then set the closed loops
            if (LoopSequences[i].IsLoop && LoopSequences[i].LType==ConcaveNonSampled)
            {
                LoopSequences[i].LType=Concave;
                Added[i]=true;
                LoopSwap.push_back(LoopSequences[i]);
            }
        }
        for (size_t i=0;i<LoopSequences.size();i++)
        {
            if (Added[i])continue;
            LoopSwap.push_back(LoopSequences[i]);
        }
        LoopSequences.clear();
        LoopSequences=LoopSwap;
    }

public:

    void GLDrawLoopInfo()
    {
        glPushAttrib(GL_ALL_ATTRIB_BITS);
        glDisable(GL_LIGHTING);
        glDepthRange(0,0.9999);


        for (size_t i=0;i<LoopSequences.size();i++)
        {
            LoopType LType= LoopSequences[i].LType;
            assert(LType!=Undetermined);

            vcg::Color4b Col=vcg::Color4b::Scatter(LoopSequences.size(),i);
            if (LType==Concave)
                Col[2]=255;
            if (LType==ConcaveNonSampled)
                Col=vcg::Color4b::Red;
            if (LType==ConvexFeature)
                Col[1]=255;

            //if (LoopSequences[i].IsLoop)
            glLineWidth(10);
            //            else
            //                glLineWidth(5);

            vcg::glColor(Col);

            glBegin(GL_LINES);
            for (size_t j=0;j<LoopSequences[i].IndexFE.size();j++)
            {
                size_t IndexF=LoopSequences[i].IndexFE[j].first;
                size_t IndexE=LoopSequences[i].IndexFE[j].second;
                assert(IndexF>=0);
                assert(IndexF<anigraph.Mesh().face.size());
                assert(IndexE>=0);
                assert(IndexE<3);
                CoordType P0=anigraph.Mesh().face[IndexF].P0(IndexE);
                CoordType P1=anigraph.Mesh().face[IndexF].P1(IndexE);
                if (LoopSequences[i].Sharp[j])
                {
                    CoordType P2=anigraph.Mesh().face[IndexF].P2(IndexE);
                    P0=P0*0.8+P2*0.2;
                    P1=P1*0.8+P2*0.2;
                }
                vcg::glVertex(P0);
                vcg::glVertex(P1);
            }
            glEnd();
        }

        glPopAttrib();
    }

    void SelectPartitionBorders(bool OnlyCrease=false)
    {
        vcg::tri::UpdateFlags<MeshType>::FaceClearFaceEdgeS(anigraph.Mesh());

        //first add borders
        for (size_t i=0;i<anigraph.Mesh().face.size();i++)
            for (size_t j=0;j<3;j++)
            {
                if (anigraph.Mesh().face[i].IsB(j))
                    anigraph.Mesh().face[i].SetFaceEdgeS(j);
            }

        for (size_t i=0;i<LoopSequences.size();i++)
            for (size_t j=0;j<LoopSequences[i].IndexFE.size();j++)
            {
                size_t IndexF=LoopSequences[i].IndexFE[j].first;
                size_t IndexE=LoopSequences[i].IndexFE[j].second;
                if (anigraph.Mesh().face[IndexF].IsB(IndexE))continue;
                if (OnlyCrease && (!LoopSequences[i].Sharp[j]))continue;

                anigraph.Mesh().face[IndexF].SetFaceEdgeS(IndexE);

                FaceType *FOpp=anigraph.Mesh().face[IndexF].FFp(IndexE);
                size_t EOpp=anigraph.Mesh().face[IndexF].FFi(IndexE);

                FOpp->SetFaceEdgeS(EOpp);
            }
    }

    void RetrievePartitioningFrom(const size_t &IndexF,std::vector<size_t> &Partition)
    {
        Partition.clear();

        std::vector<size_t> stack;
        stack.push_back(IndexF);
        do
        {
            size_t currF=stack.back();
            stack.pop_back();

            Partition.push_back(currF);
            anigraph.Mesh().face[currF].SetV();
            for (size_t j=0;j<3;j++)
            {
                if (anigraph.Mesh().face[currF].IsFaceEdgeS(j))continue;

                int NextFIndex=vcg::tri::Index(anigraph.Mesh(),anigraph.Mesh().face[currF].FFp(j));

                if (anigraph.Mesh().face[NextFIndex].IsV())continue;

                stack.push_back(NextFIndex);
            }
        }while (!stack.empty());
    }

    void GetFacePartitions(std::vector<std::vector<size_t> > &Partitions)
    {
        SelectPartitionBorders();
        Partitions.clear();
        vcg::tri::UpdateFlags<MeshType>::FaceClearV(anigraph.Mesh());
        for (size_t i=0;i<anigraph.Mesh().face.size();i++)
        {
            FaceType *f=&anigraph.Mesh().face[i];
            if (f->IsV())continue;
            Partitions.push_back(std::vector<size_t>());
            RetrievePartitioningFrom(i,Partitions.back());
        }
    }

    void SelectPathCorners()
    {
        std::vector<size_t> VertValence(anigraph.Mesh().vert.size(),0);
        std::set<std::pair<size_t,size_t> > AddedEdges;
        for (size_t i=0;i<LoopSequences.size();i++)
            for (size_t j=0;j<LoopSequences[i].IndexFE.size();j++)
            {
                size_t IndexF=LoopSequences[i].IndexFE[j].first;
                size_t IndexE=LoopSequences[i].IndexFE[j].second;
                size_t IndexV0=vcg::tri::Index(anigraph.Mesh(),anigraph.Mesh().face[IndexF].V0(IndexE));
                size_t IndexV1=vcg::tri::Index(anigraph.Mesh(),anigraph.Mesh().face[IndexF].V1(IndexE));
                std::pair<size_t,size_t> key(std::min(IndexV0,IndexV1),std::max(IndexV0,IndexV1));
                if (AddedEdges.count(key)>0)continue;
                AddedEdges.insert(key);
                VertValence[IndexV0]++;
                VertValence[IndexV1]++;
            }
        vcg::tri::UpdateFlags<MeshType>::VertexClearS(anigraph.Mesh());
        for (size_t i=0;i<anigraph.Mesh().vert.size();i++)
        {
            if (VertValence[i]==0)continue;
            if (VertValence[i]!=2)anigraph.Mesh().vert[i].SetS();
        }
    }

    void GetPartitionLoops(const std::vector<std::vector<size_t> > &Partitions,
                           std::vector<std::vector<size_t> > &PartitionsLoops)
    {
        PartitionsLoops.clear();
        PartitionsLoops.resize(Partitions.size());

        //first put all the loops in a set
        std::map<std::pair<size_t,size_t>, std::vector<size_t> > EdgeLoops;
        for (size_t i=0;i<LoopSequences.size();i++)
            for (size_t j=0;j<LoopSequences[i].IndexFE.size();j++)
            {
                size_t IndexF=LoopSequences[i].IndexFE[j].first;
                size_t IndexE=LoopSequences[i].IndexFE[j].second;
                size_t IndexV0=vcg::tri::Index(anigraph.Mesh(),anigraph.Mesh().face[IndexF].V0(IndexE));
                size_t IndexV1=vcg::tri::Index(anigraph.Mesh(),anigraph.Mesh().face[IndexF].V1(IndexE));
                std::pair<size_t,size_t> key(std::min(IndexV0,IndexV1),std::max(IndexV0,IndexV1));
                EdgeLoops[key].push_back(i);
            }

        for (size_t i=0;i<Partitions.size();i++)
        {
            for (size_t j=0;j<Partitions[i].size();j++)
            {
                size_t IndexF=Partitions[i][j];
                for (size_t j=0;j<3;j++)
                {
                    size_t IndexV0=vcg::tri::Index(anigraph.Mesh(),anigraph.Mesh().face[IndexF].V0(j));
                    size_t IndexV1=vcg::tri::Index(anigraph.Mesh(),anigraph.Mesh().face[IndexF].V1(j));
                    std::pair<size_t,size_t> key(std::min(IndexV0,IndexV1),std::max(IndexV0,IndexV1));
                    if (EdgeLoops.count(key)>0)
                        PartitionsLoops[i].insert(PartitionsLoops[i].end(),EdgeLoops[key].begin(),EdgeLoops[key].end());
                }
            }
            std::sort(PartitionsLoops[i].begin(),PartitionsLoops[i].end());
            std::vector<size_t>::iterator iteUniq=std::unique(PartitionsLoops[i].begin(),
                                                              PartitionsLoops[i].end());
            PartitionsLoops[i].erase(iteUniq, PartitionsLoops[i].end());
        }
    }

    bool OpenPartition(std::vector<size_t> &Faces,
                       std::set<std::pair<size_t,size_t> > &MustRemovedEdges)

    {
       //check the border edges
       vcg::tri::UpdateFlags<MeshType>::FaceClearV(anigraph.Mesh());
       for (size_t i=0;i<Faces.size();i++)
           anigraph.Mesh().face[Faces[i]].SetV();

       for (size_t i=0;i<Faces.size();i++)
           for (size_t j=0;j<3;j++)
           {
               size_t IndexF0=Faces[i];
               assert(anigraph.Mesh().face[IndexF0].IsV());
               size_t Edge0=j;
               size_t IndexF1=vcg::tri::Index(anigraph.Mesh(),anigraph.Mesh().face[IndexF0].FFp(Edge0));
               size_t Edge1=anigraph.Mesh().face[IndexF0].FFi(Edge0);
               //check if border
               if (anigraph.Mesh().face[IndexF1].IsV())continue;
               //check if sharp
               if (anigraph.Mesh().face[IndexF0].IsFaceEdgeS(Edge0))continue;
               assert(!anigraph.Mesh().face[IndexF1].IsFaceEdgeS(Edge1));
               MustRemovedEdges.insert(std::pair<size_t,size_t>(IndexF0,Edge0));
               MustRemovedEdges.insert(std::pair<size_t,size_t>(IndexF1,Edge1));
           }
    }

//    void GerEraseableEdges(std::)
//    {

//    }

    void RemoveLoopEdges(const std::set<std::pair<size_t,size_t> > &MustRemovedEdges)
    {
        //then remove it
        for (size_t i=0;i<LoopSequences.size();i++)
        {
            bool removed=false;
            std::vector<FaceEdgePair> IndexFESwap;
            std::vector<bool> SharpSwap;
            for (size_t j=0;j<LoopSequences[i].IndexFE.size();j++)
            {
                size_t IndexF=LoopSequences[i].IndexFE[j].first;
                size_t IndexE=LoopSequences[i].IndexFE[j].second;
                if (MustRemovedEdges.count(std::pair<size_t,size_t>(IndexF,IndexE))>0)
                {
                    removed=true;
                }
                else
                {
                    IndexFESwap.push_back(LoopSequences[i].IndexFE[j]);
                    SharpSwap.push_back(LoopSequences[i].Sharp[j]);
                }
            }
            if (removed)
            {
                LoopSequences[i].IndexFE=IndexFESwap;
                LoopSequences[i].Sharp=SharpSwap;
                LoopSequences[i].IsLoop=false;
            }
        }
    }

    void DilateLoopStep(size_t IndexLoop)
    {
        std::set<std::pair<size_t,size_t> > MustRemovedEdges;
        vcg::tri::UpdateQuality<MeshType>::VertexConstant(anigraph.Mesh(),0);
        for (size_t i=0;i<LoopSequences[IndexLoop].IndexFE.size();i++)
        {
            size_t IndexF=LoopSequences[IndexLoop].IndexFE[i].first;
            size_t IndexE=LoopSequences[IndexLoop].IndexFE[i].second;
            VertexType *v0=anigraph.Mesh().face[IndexF].V0(IndexE);
            VertexType *v1=anigraph.Mesh().face[IndexF].V1(IndexE);
            v0->Q()+=1;
            v1->Q()+=1;
        }

        for (size_t i=0;i<LoopSequences[IndexLoop].IndexFE.size();i++)
        {
            size_t IndexF=LoopSequences[IndexLoop].IndexFE[i].first;
            size_t IndexE=LoopSequences[IndexLoop].IndexFE[i].second;
            if (anigraph.Mesh().face[IndexF].IsFaceEdgeS(IndexE))continue;
            VertexType *v0=anigraph.Mesh().face[IndexF].V0(IndexE);
            VertexType *v1=anigraph.Mesh().face[IndexF].V1(IndexE);
            if ((v0->Q()==1)||(v1->Q()==1))
                MustRemovedEdges.insert(std::pair<size_t,size_t>(IndexF,IndexE));
        }
        RemoveLoopEdges(MustRemovedEdges);
    }

    void DilateLoopsStep()
    {
        SelectPartitionBorders(true);
        for (size_t i=0;i<LoopSequences.size();i++)
            if (!LoopSequences[i].IsLoop)
                DilateLoopStep(i);
    }

    void CheckPartitions()
    {
        std::vector<std::vector<size_t> > Partitions;
        GetFacePartitions(Partitions);
        std::vector<std::vector<size_t> > PartitionsLoops;
        GetPartitionLoops(Partitions,PartitionsLoops);

        std::cout<<"There are "<<Partitions.size()<<" partitions"<<std::endl;
        SelectPathCorners();

        //then select only sharp ones
        SelectPartitionBorders(true);

        //std::vector<size_t> MustSolved;
        std::set<std::pair<size_t,size_t> > MustRemovedEdges;
        for (size_t i=0;i<Partitions.size();i++)
        {
            std::set<size_t> Corners;
            for (size_t j=0;j<Partitions[i].size();j++)
            {
                size_t IndexF=Partitions[i][j];
                for (size_t j=0;j<3;j++)
                {
                    size_t IndexV=vcg::tri::Index(anigraph.Mesh(),anigraph.Mesh().face[IndexF].V(j));
                    if (!anigraph.Mesh().vert[IndexV].IsS())continue;
                    Corners.insert(IndexV);
                }
            }
            if (Corners.size()==2)
            {
                std::cout<<"WARNING: bad crossing MUST BE SOLVED "<<std::endl;
                OpenPartition(Partitions[i],MustRemovedEdges);
            }
        }
        RemoveLoopEdges(MustRemovedEdges);
        DilateLoopsStep();
        DilateLoopsStep();
    }

    void Split(const std::vector<PathType> &ChoosenPaths,
               const std::vector<bool> &IsLoop,
               const std::vector<bool> &CrossOK)
    {
        InitSideSequences();

        InitFacePathForSequence(ChoosenPaths,IsLoop);

        SplitMesh(ChoosenPaths,IsLoop);

        anigraph.Mesh().UpdateAttributes();

        //        std::cout<<"Updated Attribites"<<std::endl;
        //        vcg::tri::io::ExporterPLY<MeshType>::Save(anigraph.Mesh(),"./testM1.ply");

        //exit(0);

        std::cout<<"Updating Path"<<std::endl;
        UpdatePathForSequenceAfterSplit();

        //then color
        //ColorByPathIndex(ChoosenPaths);

        //finally retrieve the paths
        std::cout<<"Retrieving Sequences"<<std::endl;
        InitLoopSequences(CrossOK);

        std::cout<<"Done"<<std::endl;
        vcg::tri::UpdateColor<MeshType>::PerFaceConstant(anigraph.Mesh(),vcg::Color4b::White);

        CheckPartitions();
    }



    void SaveLoopInfo(std::string &PathL)
    {
        FILE *f=fopen(PathL.c_str(),"wt");
        fprintf(f,"%d\n",(int)LoopSequences.size());
        for (size_t i=0;i<LoopSequences.size();i++)
        {
            //fprintf(f,"FEATURE %d\n",(int)i);
            LoopType LType= LoopSequences[i].LType;
            assert(LType!=Undetermined);

            if (LType==Concave)
                fprintf(f,"CONCAVE\n");
            if (LType==ConcaveNonSampled)
            {
                fprintf(f,"CONCAVE\n");
                assert(!LoopSequences[i].IsLoop);
            }
            if (LType==ConvexFeature)
                fprintf(f,"CONVEX\n");
            if (LType==Regular)
                fprintf(f,"REGULAR\n");

            if (LoopSequences[i].IsLoop)
                fprintf(f,"Closed\n");
            else
                fprintf(f,"Open\n");

            if (LoopSequences[i].CrossProperly)
                fprintf(f,"Cross OK\n");
            else
                fprintf(f,"Cross FAIL\n");

            //fprintf(f,"Segment N:%d\n",(int)LoopSequences[i].IndexFE.size());
            fprintf(f,"%d\n",(int)LoopSequences[i].IndexFE.size());
            for (size_t j=0;j<LoopSequences[i].IndexFE.size();j++)
            {
                size_t IndexF=LoopSequences[i].IndexFE[j].first;
                size_t IndexE=LoopSequences[i].IndexFE[j].second;
                if (LoopSequences[i].Sharp[j])
                    fprintf(f,"%d %d 1\n",(int)IndexF,(int)IndexE);
                else
                    fprintf(f,"%d %d 0\n",(int)IndexF,(int)IndexE);
            }
        }
        fclose(f);
    }

    void GetAllFaceEdge(std::vector<std::vector<FaceEdgePair> > &FaceEdgePath)
    {
        FaceEdgePath.clear();
        FaceEdgePath.resize(LoopSequences.size());
        for (size_t i=0;i<FaceEdgePath.size();i++)
            FaceEdgePath[i]=LoopSequences[i].IndexFE;
    }

    LoopSplitter(AnisotropicGraph<MeshType> &_anigraph,
                 SharpFeaturesManager<MeshType> &_sharp) : anigraph(_anigraph),sharp(_sharp)
    {
    }
};

#endif
