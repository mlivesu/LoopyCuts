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

#ifndef SHARP_FEATURE_SAMPLER
#define SHARP_FEATURE_SAMPLER

#include <tracing_field/graph/anisotropic_graph.h>
#include <wrap/io_trimesh/export_ply.h>
#include <tracing_field/anisotropic_geodesic.h>
#include <tracing_field/loop_common_functions.h>
#include <tracing_field/sharp_feature_manager.h>
#include <vcg/complex/algorithms/mesh_to_matrix.h>

#ifndef NO_TRACING_OPENGL
#include <wrap/gl/gl_field.h>
#endif

//#include <igl/principal_curvature.h>

#define InsideW 0.0001

template < class MeshType >
class SharpFeatureSampler {
public:

    SharpFeatureSampler(AnisotropicGraph<MeshType> &_anigraph,
                        ConflictTable<MeshType> &_CTable,
                        SharpFeaturesManager<MeshType> &_sharp) : anigraph(_anigraph),CTable(_CTable),sharp(_sharp)
    {AvoidCrossIntersections=true;}

    ~SharpFeatureSampler() {}

private:
    typedef typename MeshType::CoordType    CoordType;
    typedef typename MeshType::VertexType   VertexType;
    typedef typename MeshType::FaceType     FaceType;
    typedef typename MeshType::ScalarType   ScalarType;

    AnisotropicGraph<MeshType> &anigraph;               //the mesh on which the separatrix graph is calculated
    ConflictTable<MeshType> &CTable;
    SharpFeaturesManager<MeshType> &sharp;

    typedef typename AnisotropicGraph<MeshType>::Path PathType;
    typedef typename AnisotropicGraph<MeshType>::NodeListIterator NodeListIterator;
    typedef typename AnisotropicGraph<MeshType>::NodeListConstIterator NodeListConstIterator;

    //the extreme vertices of the convex?concave patches
    //std::vector<std::vector<size_t> > ConcaveExtremes;

    //set of concave vertex for each concave feature
    std::vector<std::set<size_t> > ConcaveVert;

    //the concave features that are already closed

    //the sharp features that identifies a concave edge that need to be traced
    std::vector<std::vector<vcg::face::Pos<FaceType> > > ConcavePos;

    //those are the already solved, can be either convex or concave closed loop
    std::vector<std::vector<vcg::face::Pos<FaceType> > > GuardPos;

    //this is true if the guard is a closed concave loop or not
    std::vector<bool> IsGuardConcave;

    //the set of faces where to get the nodes to reduce the weight
    std::vector<std::vector<size_t> > Face0,Face1;

    std::vector<std::vector<size_t> > Face0Guard,Face1Guard;

    //the direction of each concave feature wrt each face
    //Note: there MUST be only one feature per face, otherwise split it
    std::vector<int> PreferredDirF;

    std::vector<int> PreferredDirGuard;

    //for each face the concave feature and the side
    //std::vector<std::pair<int,int> > PathSideF;

    //the nodes used to reduce the link weight
    std::vector<std::vector<size_t> > Nodes0,Nodes1;

    //nodes on the sharp features to avoid cross tangentially
    std::vector<std::vector<size_t> > Nodes0Guard,Nodes1Guard;

    //if a features is covered or not
    std::vector<bool > Covered0,Covered1;

    //the lenght of each feacture
    std::vector<ScalarType > LenghtFeature;

    //the vector of weight modifiers
    //std::vector<size_t> Modifier0,Modifier1;

    //the final loops on concave features
    std::vector<PathType> ConcaveLoops;
    std::vector<PathType> ProblematicLoops;

    //the vector non closed starting loops
    std::vector<size_t> NonClosedNodes;

    //nodes that cannot be used as start point
    //std::set<int> ExcludeStart;

    //the list of one with no starting nodes
    std::vector<size_t> Unsampled;

    bool AvoidCrossIntersections;
    bool CloseConvexEndPoints;
    bool CloseOnlyWithOrtho;

    ScalarType AlignFactor;
    ScalarType MaxAngle;

    std::vector<size_t> NodeConnectivityFactor;

    ScalarType TimeTrace;
    ScalarType TimeCheck;
    ScalarType TimeInit;

    //    void ColorSharpFaces()
    //    {
    //        //then colour them
    //        vcg::tri::UpdateFlags<MeshType>::FaceClearS(anigraph.Mesh());
    //        assert(Face0.size()==Face1.size());
    //        for (size_t i=0;i<Face0.size();i++)
    //        {
    //            for (size_t j=0;j<Face0[i].size();j++)
    //            {
    //                size_t FaceI=Face0[i][j];
    //                int Index=i*2;
    //                vcg::Color4b Col=vcg::Color4b::Scatter(Face0.size()*2,Index);


    //                assert(!anigraph.Mesh().face[FaceI].IsS());

    //                anigraph.Mesh().face[FaceI].SetS();
    //                anigraph.Mesh().face[FaceI].C()=Col;
    //            }
    //            for (size_t j=0;j<Face1[i].size();j++)
    //            {
    //                size_t FaceI=Face1[i][j];
    //                int Index=i*2+1;
    //                vcg::Color4b Col=vcg::Color4b::Scatter(Face1.size()*2,Index);
    //                anigraph.Mesh().face[FaceI].C()=Col;

    //                assert(!anigraph.Mesh().face[FaceI].IsS());
    //                anigraph.Mesh().face[FaceI].SetS();
    //            }
    //        }

    //        //        vcg::tri::io::ExporterPLY<MeshType>::Save(anigraph.Mesh(),"test_concave.ply",
    //        //                                                  vcg::tri::io::Mask::IOM_VERTFLAGS|
    //        //                                                  vcg::tri::io::Mask::IOM_FACECOLOR);
    //    }


    void GetFaceDirNodes(const std::vector<std::vector<size_t> > &Faces,
                         const std::vector<int> &PreferredDir,
                         std::vector<std::vector<size_t> > &Nodes)
    {
        Nodes.clear();
        Nodes.resize(Faces.size());

        for (size_t i=0;i<Faces.size();i++)
        {
            for (size_t j=0;j<Faces[i].size();j++)
            {
                int currF=Faces[i][j];
                int currD=PreferredDir[currF];
                assert(currD!=-1);
                assert(currD>=0);
                assert(currD<2);
                std::vector<size_t> nodesI;
                anigraph.faceNodesDir(currF,currD,nodesI);
                anigraph.faceNodesDir(currF,(currD+2)%4,nodesI);

                //                //get the opposite
                //                std::vector<size_t> oppI;
                //                anigraph.OppositeNode(nodesI,oppI);

                Nodes[i].insert(Nodes[i].end(),nodesI.begin(),nodesI.end());
                // Nodes0[i].insert(Nodes0[i].end(),oppI.begin(),oppI.end());
            }
            std::sort(Nodes[i].begin(),Nodes[i].end());
            std::vector<size_t>::iterator iteConc=std::unique(Nodes[i].begin(),Nodes[i].end());
            Nodes[i].erase(iteConc, Nodes[i].end());
        }
    }

public:

    void SelectConcaveSharpEndPoints()
    {
        vcg::tri::UpdateFlags<MeshType>::VertexClearS(anigraph.Mesh());
        //        for (size_t i=0;i<ConcaveExtremes.size();i++)
        //            for (size_t j=0;j<ConcaveExtremes[i].size();j++)
        //            {
        //                int IndexV=ConcaveExtremes[i][j];
        //                anigraph.Mesh().vert[IndexV].SetS();
        //            }
        std::vector<size_t> SharpVert;
        LoopFunctions<MeshType>::GetSharpEndPoints(anigraph.Mesh(),ConcavePos,SharpVert);
        for (size_t i=0;i<SharpVert.size();i++)
            anigraph.Mesh().vert[SharpVert[i]].SetS();
    }

private:

    void SelectConcaveSharpPoints()
    {
        ConcaveVert.clear();
        vcg::tri::UpdateFlags<MeshType>::VertexClearS(anigraph.Mesh());
        for (size_t i=0;i<ConcavePos.size();i++)
            for (size_t j=0;j<ConcavePos[i].size();j++)
            {
                int IndexV0=vcg::tri::Index(anigraph.Mesh(),ConcavePos[i][j].V());
                int IndexV1=vcg::tri::Index(anigraph.Mesh(),ConcavePos[i][j].VFlip());

                anigraph.Mesh().vert[IndexV0].SetS();
                anigraph.Mesh().vert[IndexV1].SetS();
            }
    }



    void InitSharpFaces()
    {
        LoopFunctions<MeshType>::GetSharpFacesInfo(anigraph.Mesh(),ConcavePos,Face0,Face1,PreferredDirF);
        GetFaceDirNodes(Face0,PreferredDirF,Nodes0);
        GetFaceDirNodes(Face1,PreferredDirF,Nodes1);

        LoopFunctions<MeshType>::GetSharpFacesInfo(anigraph.Mesh(),GuardPos,Face0Guard,Face1Guard,PreferredDirGuard);
        GetFaceDirNodes(Face0Guard,PreferredDirGuard,Nodes0Guard);
        GetFaceDirNodes(Face1Guard,PreferredDirGuard,Nodes1Guard);


        //ColorSharpFaces();
    }


    void ReduceWeightLinks(std::vector<size_t> &IndexNodes)
    {

        std::set<size_t> NodeSet(IndexNodes.begin(),IndexNodes.end());
        for (size_t i=0;i<anigraph.Graph.size();i++)
        {
            for (typename AnisotropicGraph<MeshType>::NeighborsIterator NeighIte=anigraph.Graph[i]->neighbors.begin();
                 NeighIte!=anigraph.Graph[i]->neighbors.end();NeighIte++)
            {
                int IndexN1=NeighIte->nodeID;
                if (NodeSet.count(IndexN1)==0)continue;

                //then get the index of neighbours
                (*NeighIte).Weight=InsideW;
            }
        }
    }

    void InvalidateLinks(std::vector<size_t> &IndexNodes)
    {

        std::set<size_t> NodeSet(IndexNodes.begin(),IndexNodes.end());
        for (size_t i=0;i<anigraph.Graph.size();i++)
        {
            if (NodeSet.count(i)==0)continue;
            for (typename AnisotropicGraph<MeshType>::NeighborsIterator NeighIte=anigraph.Graph[i]->neighbors.begin();
                 NeighIte!=anigraph.Graph[i]->neighbors.end();NeighIte++)
            {
                int IndexN1=NeighIte->nodeID;
                if (NodeSet.count(IndexN1)==0)continue;

                //then get the index of neighbours
                (*NeighIte).Valid=false;
            }
        }
    }

    void ReduceWeightSharpFeature()
    {
        std::vector<size_t> AllNodes;
        for (size_t i=0;i<Nodes0.size();i++)
            for (size_t j=0;j<Nodes0[i].size();j++)
                AllNodes.push_back(Nodes0[i][j]);

        for (size_t i=0;i<Nodes1.size();i++)
            for (size_t j=0;j<Nodes1[i].size();j++)
                AllNodes.push_back(Nodes1[i][j]);

        ReduceWeightLinks(AllNodes);
    }

    void InvalidateSharpFeature()
    {

        std::vector<size_t> AllNodes;
        for (size_t i=0;i<Nodes0Guard.size();i++)
            for (size_t j=0;j<Nodes0Guard[i].size();j++)
                AllNodes.push_back(Nodes0Guard[i][j]);

        for (size_t i=0;i<Nodes1Guard.size();i++)
            for (size_t j=0;j<Nodes1Guard[i].size();j++)
                AllNodes.push_back(Nodes1Guard[i][j]);

        InvalidateLinks(AllNodes);
    }

    std::vector<std::pair<CoordType,CoordType> > ReducedWLinks;

    void SetReducedConnections()
    {
        for (size_t i=0;i<anigraph.NumNodes();i++)
        {
            CoordType pos0=anigraph.Graph[i]->pos;
            typename AnisotropicGraph<MeshType>::NeighborsIterator nIt;
            for (nIt = anigraph.Graph[i]->neighbors.begin();
                 nIt != anigraph.Graph[i]->neighbors.end(); ++nIt)
            {
                if (!(*nIt).Valid)continue;
                if ((*nIt).Weight!=InsideW)continue;

                //get the neighbors node
                size_t Index1=nIt->nodeID;
                CoordType pos1=anigraph.Graph[Index1]->pos;

                ReducedWLinks.push_back(std::pair<CoordType,CoordType>(pos0,pos1));
            }
        }
    }


    bool HasSelectedVertFace(size_t IndexNode)
    {
        //get the index of the face
        size_t IndexF,M4Ind;
        anigraph.GetFaceDir(IndexNode,IndexF,M4Ind);

        //then check if the vertex is a terminal vertex of the current feature
        //it must be selected
        for (size_t i=0;i<3;i++)
            if (anigraph.Mesh().face[IndexF].V(i)->IsS())return true;
        return false;
    }


    //    void DisableFarthestEndPointNodes()
    //    {
    //        anigraph.Mesh().SelectSharpConcaveVert();
    //        //std::vector<size_t> DisableNodes;

    //        std::vector<std::vector<size_t> > testNodes;
    //        std::vector<CoordType> testPos;
    //        std::vector<ScalarType> LenghtLimit;

    //        for (size_t IndexF=0;IndexF<anigraph.Mesh().face.size();IndexF++)
    //            for (size_t IndexE=0;IndexE<3;IndexE++)
    //            {
    //                ScalarType EdgeL=(anigraph.Mesh().face[IndexF].P0(IndexE)-
    //                                  anigraph.Mesh().face[IndexF].P1(IndexE)).Norm();
    //                EdgeL/=2;
    //                EdgeL*=1.1;
    //                if (anigraph.Mesh().face[IndexF].V0(IndexE)->IsS())
    //                {
    //                    CoordType VertPos=anigraph.Mesh().face[IndexF].P0(IndexE);
    //                    for (size_t M4Ind=0;M4Ind<4;M4Ind++)
    //                    {
    //                        testPos.push_back(VertPos);
    //                        std::vector<size_t> nodesI;
    //                        anigraph.faceEdgeNodesDir(IndexF,IndexE,M4Ind,nodesI);
    //                        testNodes.push_back(nodesI);
    //                        LenghtLimit.push_back(EdgeL);
    //                    }
    //                }
    //                if (anigraph.Mesh().face[IndexF].V1(IndexE)->IsS())
    //                {
    //                    CoordType VertPos=anigraph.Mesh().face[IndexF].P1(IndexE);
    //                    for (size_t M4Ind=0;M4Ind<4;M4Ind++)
    //                    {
    //                        testPos.push_back(VertPos);
    //                        std::vector<size_t> nodesI;
    //                        anigraph.faceEdgeNodesDir(IndexF,IndexE,M4Ind,nodesI);
    //                        testNodes.push_back(nodesI);
    //                        LenghtLimit.push_back(EdgeL);
    //                    }
    //                }
    //            }
    //        assert(testNodes.size()==testPos.size());

    //        for (size_t i=0;i<testPos.size();i++)
    //        {
    //            if (testNodes[i].size()==0)continue;
    //            CoordType VertP=testPos[i];
    //            //ScalarType closeDist=std::numeric_limits<ScalarType>::max();
    //            ScalarType farthDist=0;
    //            int IndexClose=-1;
    //            for (size_t j=0;j<testNodes[i].size();j++)
    //            {
    //              ScalarType TestDist=(VertP-anigraph.NodePos(testNodes[i][j])).Norm();
    //              std::cout<<"TEst Dist:"<<TestDist<<std::endl;
    //              if (TestDist<farthDist)continue;
    //              farthDist=TestDist;
    //              IndexClose=testNodes[i][j];
    //            }
    //            assert(IndexClose!=-1);
    //            anigraph.SetValid(IndexClose,false);
    ////            for (size_t j=0;j<testNodes[i].size();j++)
    ////            {
    ////                if (testNodes[i][j]==IndexClose)continue;
    ////                anigraph.SetValid(testNodes[i][j],false);
    ////            }
    //        }
    //        anigraph.InvalidateArcsOnNonValidNodes();
    //    }

    void DisableFarthestEndPointNodes()
    {
        SelectConcaveSharpEndPoints();

        //set possible
        std::vector<int> PartitionF(anigraph.Mesh().face.size(),-1);
        for (size_t i=0;i<Face0.size();i++)
            for (size_t j=0;j<Face0[i].size();j++)
            {
                size_t IndexF=Face0[i][j];
                assert(PartitionF[IndexF]==-1);
                PartitionF[IndexF]=i;
            }

        for (size_t i=0;i<Face1.size();i++)
            for (size_t j=0;j<Face1[i].size();j++)
            {
                size_t IndexF=Face1[i][j];
                assert(PartitionF[IndexF]==-1);
                PartitionF[IndexF]=i;
            }

        std::set<std::pair<size_t,size_t>  > ExcludeFaceEdges;
        for (size_t i=0;i<ConcavePos.size();i++)
            for (size_t j=0;j<ConcavePos[i].size();j++)
            {
                vcg::face::Pos<FaceType> CurrPos=ConcavePos[i][j];
                size_t IndexF=vcg::tri::Index(anigraph.Mesh(),CurrPos.F());
                size_t IndexE=CurrPos.E();
                ExcludeFaceEdges.insert(std::pair<size_t,size_t>(IndexF,IndexE));

                CurrPos.FlipF();
                IndexF=vcg::tri::Index(anigraph.Mesh(),CurrPos.F());
                IndexE=CurrPos.E();
                ExcludeFaceEdges.insert(std::pair<size_t,size_t>(IndexF,IndexE));
            }

        std::vector<vcg::face::Pos<FaceType>  > ToFilterPos;

        //find the possible edge that must be filtered
        for (size_t IndexF=0;IndexF<anigraph.Mesh().face.size();IndexF++)
            for (size_t IndexE=0;IndexE<3;IndexE++)
            {
                std::pair<size_t,size_t> TestFaceEdge=std::pair<size_t,size_t>(IndexF,IndexE);
                if (ExcludeFaceEdges.count(TestFaceEdge)>0)continue;//on the concave itself

                //check if on partition change
                int PartitionF0=PartitionF[IndexF];
                if(PartitionF0<0)continue;
                size_t OppFIndex=vcg::tri::Index(anigraph.Mesh(),anigraph.Mesh().face[IndexF].FFp(IndexE));
                int PartitionF1=PartitionF[OppFIndex];
                if (PartitionF1==PartitionF0)continue;

                FaceType *CurrF=&anigraph.Mesh().face[IndexF];

                if ((anigraph.Mesh().face[IndexF].V0(IndexE)->IsS())
                    ||(anigraph.Mesh().face[IndexF].V1(IndexE)->IsS()))
                    ToFilterPos.push_back(vcg::face::Pos<FaceType>(CurrF,IndexE));
            }

        //push also other side
        size_t currSize=ToFilterPos.size();
        for (size_t i=0;i<currSize;i++)
        {
            vcg::face::Pos<FaceType> OppPos=ToFilterPos[i];
            OppPos.FlipF();
            ToFilterPos.push_back(OppPos);
        }

        //make them unique
        std::sort(ToFilterPos.begin(),ToFilterPos.end());
        typename std::vector<vcg::face::Pos<FaceType> >::iterator iteConc=std::unique(ToFilterPos.begin(),ToFilterPos.end());
        ToFilterPos.erase(iteConc, ToFilterPos.end());

        //put the right side (selected)
        for (size_t i=0;i<ToFilterPos.size();i++)
        {
            if (!ToFilterPos[i].V()->IsS())
                ToFilterPos[i].FlipV();
            assert(ToFilterPos[i].V()->IsS());
        }
        //        assert(testNodes.size()==testPos.size());

        for (size_t i=0;i<ToFilterPos.size();i++)
        {
            for (size_t M4Ind=0;M4Ind<4;M4Ind++)
            {
                std::vector<size_t> nodesI;
                size_t IndexF=vcg::tri::Index(anigraph.Mesh(),ToFilterPos[i].F());
                size_t IndexE=ToFilterPos[i].E();
                anigraph.faceEdgeNodesDir(IndexF,IndexE,M4Ind,nodesI);
                std::vector<std::pair<ScalarType,size_t> > NodeDist;
                CoordType TestPos=ToFilterPos[i].V()->P();

                for (size_t j=0;j<nodesI.size();j++)
                {
                    size_t TestN= nodesI[j];
                    if (!anigraph.IsValid(TestN))continue;
                    ScalarType Dist=(anigraph.NodePos(TestN)-TestPos).Norm();
                    NodeDist.push_back(std::pair<ScalarType,size_t>(Dist,TestN));
                }
                std::sort(NodeDist.begin(),NodeDist.end());

                size_t IndexN= NodeDist.back().second;
                //anigraph.SetValid(IndexN,false);
                //disable all but first
//                for (size_t j=NodeDist.size()/2;j<NodeDist.size();j++)
//                {
//                    size_t IndexN= NodeDist[j].second;
//                    //anigraph.SetValid(IndexN,false);
//                }
            }
        }
        anigraph.InvalidateArcsOnNonValidNodes();
    }

    bool MustDisabledLink(const std::vector<std::vector<int> > &NodeSet,
                          const size_t &IndexN0,const size_t &IndexN1)
    {
        bool Terminal0=HasSelectedVertFace(IndexN0);
        bool Terminal1=HasSelectedVertFace(IndexN1);
        bool IsPathNode0=(NodeSet[IndexN0].size()>0);
        bool IsPathNode1=(NodeSet[IndexN1].size()>0);

        if ((!IsPathNode0)&&(!IsPathNode1))return false;

        size_t FaceN0=anigraph.NodeFace(IndexN0);
        size_t FaceN1=anigraph.NodeFace(IndexN1);
//        bool IsInside0=(NodeSet[IndexN0].size()>0);
//        bool IsInside1=(NodeSet[IndexN1].size()>0);

        bool IsInside0=(anigraph.Mesh().face[FaceN0].IsS());
        bool IsInside1=(anigraph.Mesh().face[FaceN1].IsS());

        //both outside
        if ((!IsInside0)&&(!IsInside1))return false;

        //check terminal conditions
        if ((!IsInside0)&&(IsInside1)&&(!Terminal1))return true;
        if ((IsInside0)&&(!IsInside1)&&(!Terminal0))return true;

        for (size_t i=0;i<NodeSet[IndexN0].size();i++)
            for (size_t j=0;j<NodeSet[IndexN1].size();j++)
            {
                int Set0=NodeSet[IndexN0][i];
                int Set1=NodeSet[IndexN1][j];
                //belong to the same set
                if (Set0==Set1)continue;
                //opposite set
                if ((Set0/2)==(Set1/2))return true;
            }
        return false;
    }

    //disable when they are adjacent and not along a feature
    void DisableAdjacentFeatureFaces()
    {
        vcg::tri::UpdateSelection<MeshType>::VertexClear(anigraph.Mesh());
        LoopFunctions<MeshType>::SelectSharpEndPoints(anigraph.Mesh(),ConcavePos);

        vcg::tri::UpdateQuality<MeshType>::FaceConstant(anigraph.Mesh(),-1);
        for (size_t i=0;i<Face0.size();i++)
            for (size_t j=0;j<Face0[i].size();j++)
                anigraph.Mesh().face[Face0[i][j]].Q()=i*2;

        for (size_t i=0;i<Face1.size();i++)
            for (size_t j=0;j<Face1[i].size();j++)
                anigraph.Mesh().face[Face1[i][j]].Q()=i*2+1;

        std::vector<std::pair<size_t,size_t> > MustDisableEdges;
        std::vector<std::vector<size_t> > To_CheckFaces(Face0);
        To_CheckFaces.insert(To_CheckFaces.end(),Face1.begin(),Face1.end());
        for (size_t i=0;i<To_CheckFaces.size();i++)
            for (size_t j=0;j<To_CheckFaces[i].size();j++)
            {
                FaceType *CurrF=&anigraph.Mesh().face[To_CheckFaces[i][j]];
                int Q0=CurrF->Q();
                for (size_t e=0;e<3;e++)
                {
                    FaceType *FOpp=CurrF->FFp(e);
                    int Q1=FOpp->Q();
                    if (Q0==Q1)continue;
                    if ((Q0/2)==(Q1/2))continue;
                    if (Q0<0)continue;
                    if (Q1<0)continue;
                    if (CurrF->V0(e)->IsS())continue;
                    if (CurrF->V1(e)->IsS())continue;
                    MustDisableEdges.push_back(std::pair<size_t,size_t>(To_CheckFaces[i][j],e));
                }
            }

        for (size_t i=0;i<MustDisableEdges.size();i++)
        {
            for (size_t dir=0;dir<4;dir++)
            {
                std::vector<size_t> nodes;
                anigraph.faceEdgeNodesDir(MustDisableEdges[i].first,
                                          MustDisableEdges[i].second,
                                          dir,nodes);
                for (size_t j=0;j<nodes.size();j++)
                    anigraph.SetValid(nodes[j],false);
            }
        }
       anigraph.InvalidateArcsOnNonValidNodes();
    }

    void DisableOuterSharpLinks()
    {
        vcg::tri::UpdateSelection<MeshType>::VertexClear(anigraph.Mesh());
        LoopFunctions<MeshType>::SelectSharpEndPoints(anigraph.Mesh(),ConcavePos);
        vcg::tri::UpdateSelection<MeshType>::FaceClear(anigraph.Mesh());
        for (size_t i=0;i<Face0.size();i++)
            for (size_t j=0;j<Face0[i].size();j++)
                anigraph.Mesh().face[Face0[i][j]].SetS();
        for (size_t i=0;i<Face1.size();i++)
            for (size_t j=0;j<Face1[i].size();j++)
                anigraph.Mesh().face[Face1[i][j]].SetS();

        assert(Nodes0.size()==Nodes1.size());

        std::vector<std::vector<int> > NodeSet(anigraph.NumNodes());

        //set the partition for each node
        for (size_t i=0;i<Nodes0.size();i++)
        {
            for (size_t j=0;j<Nodes0[i].size();j++)
            {
                //assert(NodeSet[Nodes0[i][j]]==-1);
                int IndexN0=Nodes0[i][j];
                NodeSet[IndexN0].push_back(i*2);
            }
            for (size_t j=0;j<Nodes1[i].size();j++)
            {
                //assert(NodeSet[Nodes1[i][j]]==-1);
                int IndexN1=Nodes1[i][j];
                NodeSet[IndexN1].push_back((i*2)+1);
            }
        }

        for (size_t i=0;i<anigraph.NumNodes();i++)
        {
            int IndexN0=i;
            //int Set0=NodeSet[IndexN0];
            //bool Terminal0=HasSelectedVertFace(IndexN0);
            //then check all neighbours
            typename AnisotropicGraph<MeshType>::NeighborsIterator NeighIte;
            for (NeighIte=anigraph.Graph[i]->neighbors.begin();
                 NeighIte!=anigraph.Graph[i]->neighbors.end();NeighIte++)
            {
                int IndexN1=(*NeighIte).nodeID;
                bool MustDisable=MustDisabledLink(NodeSet,IndexN0,IndexN1);
                if (!MustDisable)continue;
                (*NeighIte).Valid=false;
            }
        }


        vcg::tri::UpdateSelection<MeshType>::VertexClear(anigraph.Mesh());
        vcg::tri::UpdateSelection<MeshType>::FaceClear(anigraph.Mesh());
    }


    void GetGuardFaces(const int SharpIndex,
                       const int Side,
                       std::vector<std::pair<size_t,size_t> > &FaceDirGuard)
    {
        FaceDirGuard.clear();
        if (Side==0)
        {
            for (size_t i=0;i<Face0[SharpIndex].size();i++)
            {
                int IndexF=Face0[SharpIndex][i];
                //int M4Dir=PreferredDirF[SharpIndex][IndexF];//DirF0[SharpIndex][i];
                int M4Dir=PreferredDirF[IndexF];
                FaceDirGuard.push_back(std::pair<size_t,size_t>(IndexF,M4Dir));
            }
        }else
        {
            for (size_t i=0;i<Face1[SharpIndex].size();i++)
            {
                int IndexF=Face1[SharpIndex][i];
                int M4Dir=PreferredDirF[IndexF];//DirF1[SharpIndex][i];
                FaceDirGuard.push_back(std::pair<size_t,size_t>(IndexF,M4Dir));
            }
        }
    }

    enum SampledType{completeSample,partialSample,completeUnsample};

    SampledType IsSampled(PathType CurrPath,const std::pair<int,int> &SharpSide)
    {
        int SharpIndex=SharpSide.first;
        int Side=SharpSide.second;
        assert((Side==0)||(Side==1));
        assert(SharpIndex>=0);
        assert(SharpIndex<(int)Covered0.size());
        assert(SharpIndex<(int)Covered1.size());

        //get the faces of the current Path in M2
        std::vector<std::pair<size_t,size_t> > PathFaceDir;
        anigraph.GetFacesDirM2(CurrPath,PathFaceDir,true);
        std::sort(PathFaceDir.begin(),PathFaceDir.end());
        //std::cout<<"size 0 "<<PathFaceDir.size()<<std::endl;

        //then get the set of faces and direction
        std::vector<std::pair<size_t,size_t> > FeaturesGuardFaces;
        GetGuardFaces(SharpIndex,Side,FeaturesGuardFaces);
        std::sort(FeaturesGuardFaces.begin(),FeaturesGuardFaces.end());
        //std::cout<<"size 1 "<<FeaturesGuardFaces.size()<<std::endl;


        //then check the intersection of sets
        std::vector<std::pair<size_t,size_t> > diff;
        std::set_difference(FeaturesGuardFaces.begin(), FeaturesGuardFaces.end(),
                            PathFaceDir.begin(), PathFaceDir.end(),
                            std::inserter(diff, diff.begin()));


        if (diff.empty())return completeSample;
        if (diff.size()==FeaturesGuardFaces.size())return completeUnsample;
        std::cout<<"WARNING: PARTIAL SAMPLING remaining size "<<diff.size()<<" wrt "<<FeaturesGuardFaces.size()<<std::endl;
//        ProblematicLoops.push_back(CurrPath);
        return partialSample;
    }

    struct CandidateSharpLoop
    {
        std::vector<std::pair<int,int> > CoveredFeatures;
        int StartingNode;
        std::pair<int,int> StartingPartition;
        PathType LoopP;
        ScalarType CoveredLenght;
        //ScalarType CurvatureL;
        //size_t Modifier;

        CandidateSharpLoop(const PathType &_LoopP,
                           const int _StartingNode,
                           const std::pair<int,int> _StartingPartition)
                           //ScalarType _CurvatureL)
        {
            LoopP=_LoopP;
            StartingNode=_StartingNode;
            StartingPartition=_StartingPartition;
            //CurvatureL=_CurvatureL;
            CoveredLenght=0;
            //Modifier=0;
        }

        inline bool operator <(const CandidateSharpLoop &Cloop)const
        {
            assert(StartingNode>=0);
            assert(StartingPartition.first>=0);
            assert(StartingPartition.second>=0);
            //return (CurvatureL>Cloop.CurvatureL);
            return (CoveredLenght<Cloop.CoveredLenght);

            //            int Size0=Cloop.CoveredFeatures.size();//+Cloop.Modifier;
            //            int Size1=CoveredFeatures.size();//+Modifier;

            //            if (Size0==Size1)
            //                return (LoopP.length>Cloop.LoopP.length);
            //            return (Size1<Size0);
        }
    };

    std::vector<CandidateSharpLoop> SampledLoops;

    int NumTotalUncovered()
    {
        int NumUnc=0;
        for (size_t i=0;i<Covered0.size();i++)
        {
            if (!Covered0[i])NumUnc++;
            if (!Covered1[i])NumUnc++;
        }
        return NumUnc;
    }

    //    void UpdateTotalUncoveredWeight()
    //    {
    //        for (size_t i=0;i<Covered0.size();i++)
    //        {
    //            if (!Covered0[i])Modifier0[i]+=1;
    //            if (!Covered1[i])Modifier1[i]+=1;
    //        }
    //    }

    bool IsCovered(const std::pair<int,int> &SharpSide)
    {
        int IndexSide=SharpSide.second;
        int IndexPart=SharpSide.first;
        assert((IndexSide==0)||(IndexSide==1));
        assert(IndexPart>=0);
        assert(IndexPart<(int)Covered0.size());
        assert(IndexPart<(int)Covered1.size());
        if (IndexSide==0)return (Covered0[IndexPart]);
        assert(IndexSide==1);
        return (Covered1[IndexPart]);
    }

    bool HasCoveredPartition(const CandidateSharpLoop &CLoop)
    {
        if (IsCovered(CLoop.StartingPartition))return true;
        for (size_t i=0;i<CLoop.CoveredFeatures.size();i++)
            if (IsCovered(CLoop.CoveredFeatures[i]))return true;
        return false;
    }

    bool HasCovererStartingPartition(const CandidateSharpLoop &CLoop)
    {
        if (IsCovered(CLoop.StartingPartition))return true;
        return false;
    }

    bool NeedUpdate(const CandidateSharpLoop &CLoop)
    {
        for (size_t i=0;i<CLoop.LoopP.nodes.size();i++)
        {
            int IndexN=CLoop.LoopP.nodes[i];
            if (!anigraph.IsValid(IndexN))return true;
        }
        return false;
    }


    bool UpdateCoveredFeatures(CandidateSharpLoop &CLoop)
    {
        CLoop.CoveredFeatures.clear();
        //CLoop.Modifier=0;

        //at least should cover starting partition
        SampledType STypeInitial=IsSampled(CLoop.LoopP,CLoop.StartingPartition);
        if (STypeInitial==partialSample)
        {
            std::cout<<"****WARNING: PARTIAL SAMPLE OF PARTITION ***"<<
                       CLoop.StartingPartition.first<<" , "<<
                       CLoop.StartingPartition.second<<std::endl;
            //assert(0);
            return false;
        }

        //if (STypeInitial!=completeSample){
        if (STypeInitial==completeUnsample){
            std::cout<<"WARNING: Incomplete sample of Initial Partition "<<
                       CLoop.StartingPartition.first<<" , "<<
                       CLoop.StartingPartition.second<<std::endl;
            return false;
        }

        //for each loop
        for (size_t i=0;i<Covered0.size();i++)
        {
            SampledType SType0=IsSampled(CLoop.LoopP,std::pair<int,int>(i,0));
            //if (SType0==partialSample)return false;
            SampledType SType1=IsSampled(CLoop.LoopP,std::pair<int,int>(i,1));
            //if (SType1==partialSample)return false;

            if (SType0==completeSample)
            {
                //CLoop.Modifier+=Modifier0[i];
                CLoop.CoveredFeatures.push_back(std::pair<int,int>(i,0));
                CLoop.CoveredLenght+=LenghtFeature[i];
            }
            if (SType1==completeSample)
            {
                //CLoop.Modifier+=Modifier1[i];
                CLoop.CoveredFeatures.push_back(std::pair<int,int>(i,1));
                CLoop.CoveredLenght+=LenghtFeature[i];
            }
        }
        //std::cout<<"Modifier "<<CLoop.Modifier<<std::endl;
        return true;
    }

    void UpdateTotalCoveredFeatures(CandidateSharpLoop &CLoop)
    {
        assert(CLoop.CoveredFeatures.size()>0);
        for (size_t i=0;i<CLoop.CoveredFeatures.size();i++)
        {
            int IndexSide=CLoop.CoveredFeatures[i].second;
            int IndexSharp=CLoop.CoveredFeatures[i].first;
            assert((IndexSide==0)||(IndexSide==1));
            assert(IndexSharp>=0);
            assert(IndexSharp<(int)Covered0.size());
            assert(IndexSharp<(int)Covered1.size());
            if (IndexSide==0)
            {
                assert(!Covered0[IndexSharp]);
                Covered0[IndexSharp]=true;
            }
            if (IndexSide==1)
            {
                assert(!Covered1[IndexSharp]);
                Covered1[IndexSharp]=true;
            }
        }
    }

    const size_t TargetC=5;

    void TracingNodes(const std::pair<int,int> &SharpSide,
                      std::vector<int> &StartNodes)
    {
        StartNodes.clear();
        //std::set<size_t> NodeSet;
        int IndexSide=SharpSide.second;
        int IndexSharp=SharpSide.first;
        assert((IndexSide==0)||(IndexSide==1));
        assert(IndexSharp>=0);
        assert(IndexSharp<(int)Covered0.size());
        assert(IndexSharp<(int)Covered1.size());

        std::vector<std::pair<size_t,size_t> > CandidateNodes;

        if (IndexSide==0)
        {
            for (size_t i=0;i<Nodes0[IndexSharp].size();i++)
            {
                int TestNode=Nodes0[IndexSharp][i];
                CandidateNodes.push_back(std::pair<size_t,size_t>(NodeConnectivityFactor[TestNode],TestNode));
            }
        }else
        {
            for (size_t i=0;i<Nodes1[IndexSharp].size();i++)
            {
                int TestNode=Nodes1[IndexSharp][i];
                CandidateNodes.push_back(std::pair<size_t,size_t>(NodeConnectivityFactor[TestNode],TestNode));
            }
        }
        std::sort(CandidateNodes.begin(),CandidateNodes.end());
        std::reverse(CandidateNodes.begin(),CandidateNodes.end());
        for (size_t i=0;i<std::min(TargetC,CandidateNodes.size());i++)
        {
            if (CandidateNodes[i].first==0)continue;
            StartNodes.push_back(CandidateNodes[i].second);
        }
    }


    bool AddCandidates(const std::pair<int,int> &SharpSide,std::vector<CandidateSharpLoop> &CandLoops)
    {
        std::vector<int> StartNodes;

        //get the starting tracing nodes
        TracingNodes(SharpSide,StartNodes);

        if (StartNodes.size()==0)
            Unsampled.push_back(SharpSide.first);

        int Tested=0;
        int Accepted=0;
        int NonClosed=0;
        int NonTopoOk=0;
        int CrossLoop=0;
        int NonCoverOk=0;
        bool Sampled=false;
        for (size_t i=0;i<StartNodes.size();i++)
        {
            Tested++;

            PathType CurrP;

            //std::cout<<"Tracing Loop "<<i<<std::endl;
            size_t t0=clock();
            bool closed=AnisotropicQuery<MeshType>::ExtractLoop(anigraph,StartNodes[i],MaxAngle,AlignFactor,CurrP);
            size_t t1=clock();
            TimeTrace+=(t1-t0)/(ScalarType)CLOCKS_PER_SEC;

            if (!closed)
            {
                NonClosed++;
                NonClosedNodes.push_back(StartNodes[i]);
                //std::cout<<"Warning Non Closed Nodes startin from "<<StartNodes[i]<<std::endl;
                continue;
            }
            bool SelfIntersect,CrossIntesect;
            IsTopoOkLoop(CurrP,SelfIntersect,CrossIntesect);

            if (SelfIntersect)
            {
                NonTopoOk++;
                size_t t2=clock();
                TimeCheck+=(t2-t1)/(ScalarType)CLOCKS_PER_SEC;
                //std::cout<<"Warning Topologycal OK loop startin from "<<StartNodes[i]<<std::endl;
                continue;
            }

            if ((CrossIntesect)&&(AvoidCrossIntersections))
            {
                CrossLoop++;

                size_t t2=clock();
                TimeCheck+=(t2-t1)/(ScalarType)CLOCKS_PER_SEC;

                //std::cout<<"Warning Cross intersection loop startin from "<<StartNodes[i]<<std::endl;
                continue;
            }

            //std::cout<<"Checking Loop "<<std::endl;
            //ScalarType CurrCurvLoop=LoopFunctions<MeshType>::ComputeLoopCurvature(anigraph,CurrP,PD1C,PD2C,K1,K2);
            CandidateSharpLoop CLoop(CurrP,StartNodes[i],SharpSide);//,CurrCurvLoop);

            //update the covering
            bool OkCover=UpdateCoveredFeatures(CLoop);

            //then add to candidates if is ok
            if (OkCover)
            {
                Accepted++;
            }
            else
            {
                NonCoverOk++;
                std::cout<<"WARNING: Non Covered Loops LOOP STARTING "<<StartNodes[i]<<std::endl;
                //assert(0);
                size_t t2=clock();
                TimeCheck+=(t2-t1)/(ScalarType)CLOCKS_PER_SEC;

                continue;
            }
            CandLoops.push_back(CLoop);
            Sampled=true;

            size_t t2=clock();
            TimeCheck+=(t2-t1)/(ScalarType)CLOCKS_PER_SEC;

            break;
        }
        //std::cout<<"Accepted "<<Accepted<<" out of "<<Tested<<std::endl;
        if (!Sampled)
        {
            std::cout<<"WARNING Non Sampled feature "<<std::endl;
            if (NonClosed==(int)StartNodes.size())
                std::cout<<"WARNING ALL NON CLOSED LOOPS "<<std::endl;
            if (NonTopoOk==(int)StartNodes.size())
                std::cout<<"WARNING ALL TANGENT CROSS PATHS "<<std::endl;
            if (CrossLoop==(int)StartNodes.size())
                std::cout<<"WARNING ALL ORTHO CROSS PATHS "<<std::endl;
            //            size_t t2=clock();
            //            TimeCheck+=(t2-t1)/(ScalarType)CLOCKS_PER_SEC;
        }
        return Sampled;
        //std::cout<<"Non Closed "<<NonClosed<<std::endl;
        //std::cout<<"Non Topologycal ok "<<NonTopoOk<<std::endl;
        //std::cout<<"CrossLoop "<<CrossLoop<<std::endl;
        //std::cout<<"Non Cover ok "<<NonCoverOk<<std::endl;
    }

    void SetChoosenLoopsAsBarrier()
    {
        //std::vector<bool> IsLoop(ConcaveLoops.size(),true);
        //LoopFunctions<MeshType>::BarrierPaths(anigraph,ConcaveLoops,IsLoop);

        CTable.BarrierPaths(anigraph,ConcaveLoops,false);
    }


    void InitCandidatesLoop()
    {
        LoopFunctions<MeshType>::GetNodeConnectivityFactor(anigraph,NodeConnectivityFactor);

        NonClosedNodes.clear();

        SampledLoops.clear();

        //TODO REMOVE CONVEX FEATURES!

        //initialize the covering
        Covered0=std::vector<bool>(ConcavePos.size(),false);
        Covered1=std::vector<bool>(ConcavePos.size(),false);

        size_t NumUnsampled=0;
        for (size_t i=0;i<Nodes0.size();i++)
        {
            std::cout<<"Tracing "<<i<<" side 0"<<std::endl;

            bool Sampled=AddCandidates(std::pair<int,int>(i,0),SampledLoops);

            if (!Sampled)NumUnsampled++;

            std::cout<<"Tracing "<<i<<" side 1"<<std::endl;
            Sampled=AddCandidates(std::pair<int,int>(i,1),SampledLoops);

            if (!Sampled)NumUnsampled++;
        }
        std::cout<<"THERE ARE "<<NumUnsampled<<" UNSAMPLED FEATURES"<<std::endl;
        std::cout<<"TIME "<<TimeTrace<<" FOR TRACING"<<std::endl;
        std::cout<<"TIME "<<TimeCheck<<" FOR CHECKING"<<std::endl;
    }


    void IsTopoOkLoop(PathType &P,bool &SelfConflict,bool &CrossIntersect)
    {
        SelfConflict=false;
        CrossIntersect=false;

        ConflictFinder<MeshType> ConflFinder(anigraph);
        SelfConflict=ConflFinder.SelfConflict(P,2);
        if (SelfConflict)return;

        if (AvoidCrossIntersections)
            CrossIntersect=ConflFinder.CrossIntersect(P,true);
    }

    void UpdateCandidatesLoop()
    {
        std::vector<CandidateSharpLoop> SwapCand;

        //set the choosen loop as barrier
        SetChoosenLoopsAsBarrier();

        int Retraced=0;

        std::cout<<"Updating Step"<<std::endl;

        //then re-trace the ones need to be retracen
        for (size_t i=0;i<SampledLoops.size();i++)
        {
            //in such case the sharp feature has already been sampled
            //if (!IsOkCandidate(SampledLoops[i]))continue;
            if (HasCovererStartingPartition(SampledLoops[i]))continue;
            if (NeedUpdate(SampledLoops[i]))
            {
                //retrace it
                CandidateSharpLoop Cloop=SampledLoops[i];
                //                int PartitionI=Cloop.StartingPartition.first;
                //                int SideI=Cloop.StartingPartition.second;
                Retraced++;
                AddCandidates(Cloop.StartingPartition,SwapCand);
                //                int StartNode=Cloop.StartingNode;
                //                PathType CurrP;

                //                Retraced++;

                //                //then do the actual tracing
                //                bool closed=AnisotropicQuery<MeshType>::ExtractLoop(anigraph,StartNode,45,AlignFactor,CurrP);
                //                if (!closed)continue;
                //                bool SelfIntersect,CrossIntesect;
                //                IsTopoOkLoop(CurrP,SelfIntersect,CrossIntesect);
                //                if (SelfIntersect || CrossIntesect)continue;

                //                //then update the loop
                //                Cloop.LoopP=CurrP;

                //                //then update covered features
                //                UpdateCoveredFeatures(Cloop);

                //                //the check if cover its partition
                //                if (!IsOkCandidate(Cloop))continue;

                //otherwise add the current
                //                SwapCand.push_back(Cloop);
            }
            else
                SwapCand.push_back(SampledLoops[i]);
        }
        SampledLoops.clear();
        SampledLoops=SwapCand;
        std::sort(SampledLoops.begin(),SampledLoops.end());
        std::cout<<"Retraced "<<Retraced<<" loops"<<std::endl;
    }

    void ColorSharpByCovering()
    {
        for (size_t i=0;i<Face0.size();i++)
        {
            vcg::Color4b Col0(0,255,0,255);
            if (!Covered0[i])Col0=vcg::Color4b(255,0,0,255);
            vcg::Color4b Col1(0,255,0,255);
            if (!Covered1[i])Col1=vcg::Color4b(255,0,0,255);

            for (size_t j=0;j<Face0[i].size();j++)
                anigraph.Mesh().face[Face0[i][j]].C()=Col0;

            for (size_t j=0;j<Face1[i].size();j++)
                anigraph.Mesh().face[Face1[i][j]].C()=Col1;
        }
    }

    int TraceLoops()
    {
        SampledLoops.clear();
        ConcaveLoops.clear();

        //initialize candidates
        InitCandidatesLoop();

        //        for (size_t i=0;i<SampledLoops.size();i++)
        //            ConcaveLoops.push_back(SampledLoops[i].LoopP);

        //sort them
        std::sort(SampledLoops.begin(),SampledLoops.end());

        //        ConcaveLoops.clear();

        int NumUnsampled=NumTotalUncovered();

        if (!SampledLoops.empty())
            do{
            if (HasCoveredPartition(SampledLoops.back()))
            {
                SampledLoops.pop_back();
            }else
            {
                //assert(!HasCoveredPartition(SampledLoops.back()));
                //check if other ones
                //            if (IsOkCandidate(SampledLoops.back()))
                //            {
                ConcaveLoops.push_back(SampledLoops.back().LoopP);

                UpdateTotalCoveredFeatures(SampledLoops.back());

                //remove the element on the stack
                SampledLoops.pop_back();

                std::cout<<"Added loop"<<std::endl;

                //and update loops
                UpdateCandidatesLoop();

                //            }
                //            else //otherwise juse remove it
                //                SampledLoops.pop_back();

                NumUnsampled=NumTotalUncovered();
            }
        }while ((!SampledLoops.empty())&&(NumUnsampled>0));

        //UpdateTotalUncoveredWeight();


        int ret=NumTotalUncovered();
        return ret;
    }



    void GetConcaveNodes(std::vector<size_t> &ConcaveNodes,
                         bool Opposite=false)
    {
        ConcaveNodes.clear();
        for (size_t i=0;i<Nodes0.size();i++)
            for (size_t j=0;j<Nodes0[i].size();j++)
            {
                ConcaveNodes.push_back(Nodes0[i][j]);
                if (!Opposite)continue;
                ConcaveNodes.push_back(anigraph.OppositeNode(Nodes0[i][j]));
            }

        for (size_t i=0;i<Nodes1.size();i++)
            for (size_t j=0;j<Nodes1[i].size();j++)
            {
                ConcaveNodes.push_back(Nodes1[i][j]);
                if (!Opposite)continue;
                ConcaveNodes.push_back(anigraph.OppositeNode(Nodes1[i][j]));
            }

        std::sort(ConcaveNodes.begin(),ConcaveNodes.end());
        std::vector<size_t>::iterator iteConc=std::unique(ConcaveNodes.begin(),ConcaveNodes.end());
        ConcaveNodes.erase(iteConc, ConcaveNodes.end());
    }

    void GetExcludeDisablingNode(std::vector<size_t> &ExcludeNodes)
    {
        //select concave end points
        //anigraph.Mesh().SelectEndSharpVert(ConcaveEdge);
        sharp.SelectEndSharpVert(ConcaveEdge);

        std::vector<std::vector<std::pair<int,int> > > VertFaceDir;
        VertFaceDir.resize(anigraph.Mesh().vert.size());
        for (size_t i=0;i<ConcavePos.size();i++)
            for (size_t j=0;j<ConcavePos[i].size();j++)
            {
                VertexType *v0=ConcavePos[i][j].V();
                VertexType *v1=ConcavePos[i][j].VFlip();
                assert(v0!=v1);
                int IndexF=vcg::tri::Index(anigraph.Mesh(),ConcavePos[i][j].F());
                int PrefDir=PreferredDirF[IndexF];
                assert(PrefDir!=-1);
                if (v0->IsS())
                {
                    std::pair<int,int> Info(IndexF,PrefDir);
                    int IndxV0=vcg::tri::Index(anigraph.Mesh(),v0);
                    VertFaceDir[IndxV0].push_back(Info);
                }
                if (v1->IsS())
                {
                    std::pair<int,int> Info(IndexF,PrefDir);
                    int IndxV1=vcg::tri::Index(anigraph.Mesh(),v1);
                    VertFaceDir[IndxV1].push_back(Info);
                }
            }

        //then cycle over all faces
        std::vector<std::pair<int,int> > ExcludeFaceDir;
        for (size_t i=0;i<anigraph.Mesh().face.size();i++)
        {
            for (int j=0;j<anigraph.Mesh().face[i].VN();j++)
            {
                if (!anigraph.Mesh().face[i].V(j)->IsS())continue;
                int IndexV=vcg::tri::Index(anigraph.Mesh(),anigraph.Mesh().face[i].V(j));

                //assert(VertFaceDir[IndexV].size()>0);
                for (size_t k=0;k<VertFaceDir[IndexV].size();k++)
                {
                    int IndexF=VertFaceDir[IndexV][k].first;
                    int M4Dir=VertFaceDir[IndexV][k].second;
                    //the face itself
                    if (IndexF==(int)i)continue;
                    int NewDir=vcg::tri::CrossField<MeshType>::FollowDirection(anigraph.Mesh().face[IndexF],
                                                                               anigraph.Mesh().face[i],M4Dir);
                    ExcludeFaceDir.push_back(std::pair<int,int>(i,NewDir));
                }
            }
        }

        ExcludeNodes.clear();
        for (size_t i=0;i<Nodes0.size();i++)
        {
            ExcludeNodes.insert(ExcludeNodes.end(),Nodes0[i].begin(),Nodes0[i].end());
            ExcludeNodes.insert(ExcludeNodes.end(),Nodes1[i].begin(),Nodes1[i].end());
        }

        for (size_t i=0;i<ExcludeFaceDir.size();i++)
        {
            std::vector<size_t> nodesI;
            int currF=ExcludeFaceDir[i].first;
            int currD=ExcludeFaceDir[i].second;
            currD=currD%2;
            anigraph.faceNodesDir(currF,currD,nodesI);
            anigraph.faceNodesDir(currF,currD+2,nodesI);
            ExcludeNodes.insert(ExcludeNodes.end(),nodesI.begin(),nodesI.end());
            //std::cout<<"DE"<<std::endl;
        }
    }

public:

    void InitStructures()
    {
        TimeInit=0;
        size_t t0=clock();

        ConcaveVert.clear();
        ConcavePos.clear();
        GuardPos.clear();
        LenghtFeature.clear();

        IsGuardConcave.clear();
        Face0.clear();
        Face1.clear();
        Face0Guard.clear();
        Face1Guard.clear();
        PreferredDirF.clear();
        PreferredDirGuard.clear();

        Nodes0.clear();
        Nodes1.clear();

        Nodes0Guard.clear();
        Nodes1Guard.clear();

        Covered0.clear();
        Covered1.clear();

        //        Modifier0.clear();
        //        Modifier1.clear();

        ConcaveLoops.clear();

        NonClosedNodes.clear();

        //ExcludeStart.clear();

        vcg::tri::UpdateFlags<MeshType>::FaceClear(anigraph.Mesh());
        vcg::tri::UpdateFlags<MeshType>::VertexClear(anigraph.Mesh());
        anigraph.Mesh().UpdateAttributes();

        std::vector<std::vector<vcg::face::Pos<FaceType> > > FeatureSeq;
        std::vector<SharpFeatureType> FeatureType;
        //anigraph.Mesh().GetSharpFeatures(FeatureSeq,FeatureType);
        sharp.GetSharpFeatures(FeatureSeq,FeatureType);

        for (size_t i=0;i<FeatureType.size();i++)
        {
            if (FeatureType[i]==ConvexEdge)
            {
                GuardPos.push_back(FeatureSeq[i]);
                IsGuardConcave.push_back(false);
            }
            else
            {
                //check if it is closed or not
                if (LoopFunctions<MeshType>::IsClosedLoop(FeatureSeq[i]))
                {
                    //no need to trace but keep the info
                    GuardPos.push_back(FeatureSeq[i]);
                    IsGuardConcave.push_back(true);
                }
                else
                {
                    //otherwise need to trace it
                    ConcavePos.push_back(FeatureSeq[i]);
                    LenghtFeature.push_back(LoopFunctions<MeshType>::PosLenght(FeatureSeq[i]));
                }
            }
        }

        std::cout<<"There are "<<GuardPos.size()<<" guard pos features "<<std::endl;
        std::cout<<"There are "<<ConcavePos.size()<<" concave features "<<std::endl;


        InitSharpFaces();

        ReduceWeightSharpFeature();

        anigraph.SetAllValid();

        DisableOuterSharpLinks();
        DisableAdjacentFeatureFaces();

        //TEST NEW FUNCTION
        DisableFarthestEndPointNodes();

        std::vector<size_t> ExcludeNodes;
        GetExcludeDisablingNode(ExcludeNodes);

        //        for (size_t i=0;i<Nodes0.size();i++)
        //        {
        //            ExcludeNodes.insert(ExcludeNodes.end(),Nodes0[i].begin(),Nodes0[i].end());
        //            ExcludeNodes.insert(ExcludeNodes.end(),Nodes1[i].begin(),Nodes1[i].end());
        //        }

        //InvalidateSharpFeature();

        LoopFunctions<MeshType>::SetBarrierSharp(anigraph,GuardPos,ExcludeNodes);



        SetReducedConnections();


        size_t t1=clock();
        TimeInit=(t1-t0)/(ScalarType)CLOCKS_PER_SEC;

        //InitExcludeSamplingNodes();
    }

#ifndef NO_TRACING_OPENGL
    void GLDrawFeature(std::vector<vcg::face::Pos<FaceType> > &PathPos)
    {
        glPushAttrib(GL_ALL_ATTRIB_BITS);
        glDisable(GL_LIGHTING);
        glDepthRange(0,0.99);


        glLineWidth(20);
        glBegin(GL_LINES);

        for (size_t j=0;j<PathPos.size();j++)
        {
            FaceType *f=PathPos[j].F();
            int IndexE=PathPos[j].E();
            CoordType Pos0=f->P0(IndexE);
            CoordType Pos1=f->P1(IndexE);
            vcg::glVertex(Pos0);
            vcg::glVertex(Pos1);
        }

        glEnd();
        glPopAttrib();
    }
#endif
    bool IsSelectedPosV(vcg::face::Pos<FaceType> &Pos)
    {
        //see if has a selected vertex at the beginning
        if (Pos.V()->IsS())
        {
            assert(!Pos.VFlip()->IsS());
            return true;
        }
        if (Pos.VFlip()->IsS())
        {
            assert(!Pos.V()->IsS());
            Pos.FlipV();
            return true;
        }
        return false;
    }

    void GetTargetNodesOppositeSharpEndPos(const vcg::face::Pos<FaceType> &StartPos,
                                           const std::vector<size_t> &ConnFactor,
                                           int &N0, int &N1)
    {
        vcg::face::Pos<FaceType> TestPos=StartPos;
        CoordType TestDir=TestPos.VFlip()->P()-TestPos.V()->P();
        TestDir.Normalize();

        size_t E0,E1,F0,F1;
        TestPos.FlipE();
        E0=TestPos.E();
        TestPos.FlipE();
        F0=vcg::tri::Index(anigraph.Mesh(),TestPos.F());
        TestPos.FlipF();
        F1=vcg::tri::Index(anigraph.Mesh(),TestPos.F());
        TestPos.FlipE();
        E1=TestPos.E();

        std::vector<size_t> nodesI0;
        anigraph.faceEdgeNodesDirCoord(F0,E0,TestDir,nodesI0);

        std::vector<size_t> nodesI1;
        anigraph.faceEdgeNodesDirCoord(F1,E1,TestDir,nodesI1);

        //std::cout<<"A"<<std::endl;
        N0=-1;
        N1=-1;
        size_t MaxC0=0;
        size_t MaxC1=0;
        for (size_t i=0;i<nodesI0.size();i++)
        {
            //std::cout<<"B"<<std::endl;
            if (ConnFactor[nodesI0[i]]<MaxC0)continue;
            MaxC0=ConnFactor[nodesI0[i]];
            N0=nodesI0[i];
        }

        for (size_t i=0;i<nodesI1.size();i++)
        {
            if (ConnFactor[nodesI1[i]]<MaxC1)continue;
            MaxC1=ConnFactor[nodesI1[i]];
            N1=nodesI1[i];
        }
    }

    std::vector<int> TargetN;
    std::vector<int> SourceN;
    std::vector<bool> UsedPortal;
    typedef std::pair<int,int> FeatureSideType;
    std::vector<FeatureSideType> Feature_Side_Index;
    //    std::vector<bool> Feature_Side_Sampled;
    std::vector<PathType> ConvexCompletionPaths;
    std::vector<PathType> ConvexCompletionLoops;

    std::set<std::pair<FeatureSideType,FeatureSideType> > InsertedPaths;

    //std::vector
    void InitIncompleteConvexNodes()
    {
        Feature_Side_Index.clear();

        std::vector<size_t> ConnFactor;
        LoopFunctions<MeshType>::GetNodeConnectivityFactor(anigraph,ConnFactor);

        //first get convex ones
        std::vector<std::vector<vcg::face::Pos<FaceType> > > FeatureSeq;
        std::vector<SharpFeatureType> FeatureType;
        //anigraph.Mesh().GetSharpFeatures(FeatureSeq,FeatureType);
        sharp.GetSharpFeatures(FeatureSeq,FeatureType);

        std::vector<std::vector<vcg::face::Pos<FaceType> > > ConvexPos;

        for (size_t i=0;i<FeatureType.size();i++)
        {
            if (FeatureType[i]!=ConvexEdge)continue;
            ConvexPos.push_back(FeatureSeq[i]);
        }

        //get the ones that have incomplete end points
        //size_t NumSel=anigraph.Mesh().SelectEndSharpVert();

        //size_t NumSel=anigraph.Mesh().SelectEndSharpVert();
        size_t NumSel=sharp.SelectEndSharpVert();

        //size_t NumSel=anigraph.Mesh().SelectValence4EndPos();
        if (NumSel==0)return;

        for (size_t i=0;i<ConvexPos.size();i++)
        {
            vcg::face::Pos<FaceType> TestPos=ConvexPos[i][0];
            if (IsSelectedPosV(TestPos))
            {
                int N0,N1;
                GetTargetNodesOppositeSharpEndPos(TestPos,ConnFactor,N0,N1);

                Feature_Side_Index.push_back(std::pair<int,int>(i,0));
                Feature_Side_Index.push_back(std::pair<int,int>(i,0));

                assert(N0>=0);
                assert(N1>=0);

                TargetN.push_back(N0);
                TargetN.push_back(N1);

                int N0Opp=-1;
                int N1Opp=-1;

                N0Opp=anigraph.OppositeNode(N0);
                N1Opp=anigraph.OppositeNode(N1);

                assert(N0Opp>=0);
                assert(N1Opp>=0);

                SourceN.push_back(N0Opp);
                SourceN.push_back(N1Opp);
            }
            TestPos=ConvexPos[i].back();
            if (IsSelectedPosV(TestPos))
            {
                int N0,N1;
                GetTargetNodesOppositeSharpEndPos(TestPos,ConnFactor,N0,N1);

                Feature_Side_Index.push_back(std::pair<int,int>(i,1));
                Feature_Side_Index.push_back(std::pair<int,int>(i,1));

                assert(N0>=0);
                assert(N1>=0);
                TargetN.push_back(N0);
                TargetN.push_back(N1);

                int N0Opp=-1;
                int N1Opp=-1;
                N0Opp=anigraph.OppositeNode(N0);
                N1Opp=anigraph.OppositeNode(N1);
                assert(N0Opp>=0);
                assert(N1Opp>=0);

                SourceN.push_back(N0Opp);
                SourceN.push_back(N1Opp);
            }
        }

        UsedPortal.resize(TargetN.size(),false);
        InsertedPaths.clear();
    }

    struct ConvexSampleCandidate
    {
        PathType Path;
        FeatureSideType Side0;
        FeatureSideType Side1;
        int IndexS;
        int IndexT;

        void CheckValid()const
        {
            assert(Side0.first>=0);
            assert(Side0.second>=0);
            assert(Side1.first>=0);
            assert(Side1.second>=0);
            assert(Side0!=Side1);
        }

        inline bool operator <(const ConvexSampleCandidate &c1)const
        {
            CheckValid();
            c1.CheckValid();
            return (Path.length>c1.Path.length);
        }

        ConvexSampleCandidate()
        {
            Side0=std::pair<int,int>(-1,-1);
            Side1=std::pair<int,int>(-1,-1);
        }

        ConvexSampleCandidate(PathType &_Path,
                              std::pair<int,int> &_Side0,
                              std::pair<int,int> &_Side1,
                              int &_IndexS,
                              int &_IndexT)
        {
            Path=_Path;
            Side0=std::min(_Side0,_Side1);
            Side1=std::max(_Side0,_Side1);
            IndexS=_IndexS;
            IndexT=_IndexT;
            CheckValid();
        }
    };

    std::vector<ConvexSampleCandidate> ConvexCandidates;


    bool TracePath(std::vector<size_t> IndexSources,
                   int IndexTarget,
                   PathType &Path)
    {
        assert(IndexTarget>=0);
        assert(IndexTarget<(int)anigraph.NumNodes());
        typename AnisotropicQuery<MeshType>::DijkstraParam Param;
        Param.sourceGroups.push_back(IndexSources);
        Param.target=IndexTarget;
        int Node=AnisotropicQuery<MeshType>::computePerNodeDijsktra(anigraph,Param,MaxAngle,AlignFactor);
        if (Node==-1)return false;

        AnisotropicQuery<MeshType>::retrievePath(anigraph,(size_t)IndexTarget,Path);
        return true;
    }

    bool InitConvexSampleCandidate()
    {
        //initialize the map for each source node
        std::map<int,std::pair<int,int> > FeatureNode;
        for (size_t i=0;i<TargetN.size();i++)
            FeatureNode[TargetN[i]]=Feature_Side_Index[i];
        for (size_t i=0;i<SourceN.size();i++)
            FeatureNode[SourceN[i]]=Feature_Side_Index[i];

        ConvexCandidates.clear();
        for (size_t i=0;i<TargetN.size();i++)
        {
            if (UsedPortal[i])continue;

            int IndexTarget=i;

            int CurrTarget=TargetN[i];
            std::pair<int,int> FeatureSide0=Feature_Side_Index[i];
            assert(FeatureSide0.first>=0);
            assert(FeatureSide0.second>=0);

            std::vector<size_t> IndexSources;
            for (size_t j=0;j<SourceN.size();j++)
            {
                if (UsedPortal[j])continue;

                int CurrSource=SourceN[j];
                std::pair<int,int> FeatureSideTest=Feature_Side_Index[j];
                if (FeatureSide0==FeatureSideTest)continue;
                IndexSources.push_back(CurrSource);
            }
            //std::cout<<"Tracing Convex "<<i<<std::endl;
            PathType Path;
            bool Traced=TracePath(IndexSources,CurrTarget,Path);
            if (!Traced)continue;

            //retrieve feature side
            std::pair<int,int> FeatureSide1(-1,-1);
            int IndexSource=-1;
            for (size_t j=0;j<SourceN.size();j++)
            {
                if (SourceN[j]!=(int)Path.nodes[0])continue;
                FeatureSide1=FeatureNode[SourceN[j]];
                IndexSource=j;
            }

            assert(FeatureSide1.first>=0);
            assert(FeatureSide1.second>=0);


            ConvexSampleCandidate CSample(Path,FeatureSide0,FeatureSide1,IndexSource,IndexTarget);
            std::pair<FeatureSideType,FeatureSideType> Entry(CSample.Side0,CSample.Side1);
            if (InsertedPaths.count(Entry)>0)continue;//already done

            ConvexCandidates.push_back(CSample);
        }
        std::sort(ConvexCandidates.begin(),ConvexCandidates.end());
        std::cout<<"Size Candidates: "<<ConvexCandidates.size()<<std::endl;
        return(ConvexCandidates.size()>0);
    }

    void AddSharpConvexFeature(ConvexSampleCandidate &CSample)
    {
        ConvexCompletionPaths.push_back(CSample.Path);

        //then add the sampled element
        int Index0=CSample.IndexS;
        int Index1=CSample.IndexT;

        //update the sampled one
        UsedPortal[Index0]=true;
        UsedPortal[Index1]=true;

        std::pair<FeatureSideType,FeatureSideType> Entry(CSample.Side0,CSample.Side1);
        assert(InsertedPaths.count(Entry)==0);
        InsertedPaths.insert(Entry);
        CTable.AddBarrierNonLoop(anigraph,CSample.Path);
    }

    void TraceConvexCompletionPaths()
    {

        ConvexCompletionPaths.clear();
        bool HasCand=InitConvexSampleCandidate();
        while (HasCand)
        {
            AddSharpConvexFeature(ConvexCandidates.back());
            HasCand=InitConvexSampleCandidate();
        }
    }

    size_t GetOrthoTracing(int IndexN,const std::vector<size_t> &ConnFactor)
    {
        size_t faceId,edgeId,M4Dir;
        std::vector<size_t> nodesI;

        anigraph.GetFaceEdgeDir(IndexN,faceId,edgeId,M4Dir);

        M4Dir=(M4Dir+1)%4;
        anigraph.faceEdgeNodesDir(faceId,edgeId,M4Dir,nodesI);
        M4Dir=(M4Dir+2)%4;
        anigraph.faceEdgeNodesDir(faceId,edgeId,M4Dir,nodesI);

        size_t MaxC=0;
        size_t RetNode=0;
        for (size_t i=0;i<nodesI.size();i++)
        {
            if (ConnFactor[nodesI[i]]<MaxC)continue;
            MaxC=ConnFactor[nodesI[i]];
            RetNode=nodesI[i];
        }
        return RetNode;
    }

    void TraceIncompletedEndPos()
    {
        std::vector<size_t> ConnFactor;
        LoopFunctions<MeshType>::GetNodeConnectivityFactor(anigraph,ConnFactor);
        std::set<FeatureSideType> SideFeatures;
        std::vector<size_t> StartNodeLoops;
        for (size_t i=0;i<UsedPortal.size();i++)
        {
            if (!UsedPortal[i])continue;
            SideFeatures.insert(Feature_Side_Index[i]);
        }

        for (size_t i=0;i<UsedPortal.size();i++)
        {
            if (UsedPortal[i])continue;
            if (SideFeatures.count(Feature_Side_Index[i])>0)continue;
            SideFeatures.insert(Feature_Side_Index[i]);
            int IndexN=SourceN[i];
            StartNodeLoops.push_back(GetOrthoTracing(IndexN,ConnFactor));
        }
        std::cout<<"Remaining to trace "<<StartNodeLoops.size()<<std::endl;

        ConvexCompletionLoops.clear();
        for (size_t i=0;i<StartNodeLoops.size();i++)
        {
            PathType CurrP;
            bool closed=AnisotropicQuery<MeshType>::ExtractLoop(anigraph,StartNodeLoops[i],MaxAngle,AlignFactor,CurrP);
            if (!closed)continue;
                        bool SelfIntersect,CrossIntesect;
                        IsTopoOkLoop(CurrP,SelfIntersect,CrossIntesect);
                        if (SelfIntersect)continue;
                        if ((CrossIntesect)&&(AvoidCrossIntersections))continue;
            ConvexCompletionLoops.push_back(CurrP);
        }
    }

public:

    void TestConvexClose()
    {
        InitIncompleteConvexNodes();
        TraceConvexCompletionPaths();
        TraceIncompletedEndPos();
        SnapLoops();
    }

    void ConvexClose()
    {
        InitIncompleteConvexNodes();

        if (!CloseOnlyWithOrtho)
            TraceConvexCompletionPaths();

        TraceIncompletedEndPos();
        //SnapLoops();
    }

    void InitCandidatesLoopTest(ScalarType _AlignFactor,
                                ScalarType _MaxAngle,
                                bool _AvoidCrossIntersections)
    {
        AvoidCrossIntersections=_AvoidCrossIntersections;

        AlignFactor=_AlignFactor;
        MaxAngle=_MaxAngle;

        //        Modifier0=std::vector<size_t>(Nodes0.size(),0);
        //        Modifier1=std::vector<size_t>(Nodes1.size(),0);

        SampledLoops.clear();
        NonClosedNodes.clear();

        //initialize the covering
        Covered0=std::vector<bool>(ConcavePos.size(),false);
        Covered1=std::vector<bool>(ConcavePos.size(),false);


        TimeTrace=0;
        TimeCheck=0;
        InitCandidatesLoop();

    }

    void SampleLoopTest(ScalarType _AlignFactor,bool _AvoidCrossIntersections)
    {
        AvoidCrossIntersections=_AvoidCrossIntersections;

        AlignFactor=_AlignFactor;

        AlignFactor=_AlignFactor;
        AvoidCrossIntersections=_AvoidCrossIntersections;

        //LoopFunctions<MeshType>::CheckOppositeValid(anigraph);
        InitStructures();

        TimeTrace=0;
        TimeCheck=0;

        TraceLoops();

    }

#ifndef NO_TRACING_OPENGL
    void GLDrawSharpDebug()
    {
        ScalarType sizeDiag=0.0002;
        for (size_t i=0;i<TargetN.size();i++)
        {
            int IndexN=TargetN[i];
            if (IndexN>=0)
                anigraph.GlDrawNodeDebugInfo(IndexN,MaxAngle,AlignFactor,sizeDiag,true,false,vcg::Color4b::Red);
        }
        for (size_t i=0;i<SourceN.size();i++)
        {
            int IndexN=SourceN[i];
            if (IndexN>=0)
                anigraph.GlDrawNodeDebugInfo(IndexN,MaxAngle,AlignFactor,sizeDiag,true,false,vcg::Color4b::Blue);
        }
        //        for (size_t i=0;i<ConvexCandidates.size();i++)
        //        {
        //            anigraph.GlDrawPath(ConvexCandidates[i].Path,vcg::Color4b::Green,5,false);
        //        }
        for (size_t i=0;i<ConvexCompletionPaths.size();i++)
            anigraph.GlDrawPath(ConvexCompletionPaths[i],vcg::Color4b::Green,20,false);

        for (size_t i=0;i<ConvexCompletionLoops.size();i++)
            anigraph.GlDrawPath(ConvexCompletionLoops[i],vcg::Color4b::Green,20,true);

    }

    void GLDrawDebug(bool DrawNode,
                     bool DrawNeigh,
                     bool DrawDisabled,
                     bool DrawCandidates)
    {
        ScalarType sizeDiag=0.0001;

        glPushAttrib(GL_ALL_ATTRIB_BITS);
        glDisable(GL_LIGHTING);
        glDepthRange(0,0.99999);
        vcg::glColor(vcg::Color4b(255,0,255,255));
        glLineWidth(10);
        glBegin(GL_LINES);
        for (size_t i=0;i<PreferredDirF.size();i++)
        {
            int IndexDir=PreferredDirF[i];
            if (IndexDir==-1)continue;
            CoordType P0=(anigraph.Mesh().face[i].P(0)+
                          anigraph.Mesh().face[i].P(1)+
                          anigraph.Mesh().face[i].P(2))/3;
            CoordType DirF;
            if (IndexDir==0)
                DirF=anigraph.Mesh().face[i].PD1();
            else
                DirF=anigraph.Mesh().face[i].PD2();
            //std::cout<<"A"<<std::endl;
            DirF*=anigraph.Mesh().bbox.Diag()*sizeDiag*10;
            CoordType P1=P0+DirF;
            vcg::glVertex(P0);
            vcg::glVertex(P1);
        }
        glEnd();
        glPopAttrib();

        if (DrawNode)
        {
            vcg::glColor(vcg::Color4b(255,0,0,255));
            for (size_t i=0;i<Unsampled.size();i++)
            {
                //std::cout<<"DEEEH BOIA"<<std::endl;
                //std::cout<<ConcavePos[Unsampled[i]].size()<<std::endl;
                GLDrawFeature(ConcavePos[Unsampled[i]]);
            }

            for (size_t i=0;i<NonClosedNodes.size();i++)
            {
                //std::cout<<"DEEEH"<<std::endl;
                size_t IndexN=NonClosedNodes[i];
                anigraph.GlDrawNodeDebugInfo(IndexN,MaxAngle,AlignFactor,sizeDiag,true,true,vcg::Color4b::Red);
            }
            //            for (size_t i=0;i<Nodes0.size();i++)
            //                for (size_t j=0;j<Nodes0[i].size();j++)
            //                    anigraph.GlDrawNodeDebugInfo(Nodes0[i][j],45,20,sizeDiag,true,true,vcg::Color4b::Green);

            //            for (size_t i=0;i<Nodes1.size();i++)
            //                for (size_t j=0;j<Nodes1[i].size();j++)
            //                    anigraph.GlDrawNodeDebugInfo(Nodes1[i][j],45,20,sizeDiag,true,true,vcg::Color4b::Green);

        }



        if (DrawNeigh)
        {
            glPushAttrib(GL_ALL_ATTRIB_BITS);
            glDisable(GL_LIGHTING);
            glDepthRange(0,0.99999);

            vcg::glColor(vcg::Color4b(0,255,0,255));
            glLineWidth(1);
            glBegin(GL_LINES);

            for (size_t i=0;i<ReducedWLinks.size();i++)
            {
                vcg::glVertex(ReducedWLinks[i].first);
                vcg::glVertex(ReducedWLinks[i].second);
            }

            glEnd();
            glPopAttrib();
        }
        if (DrawDisabled)
        {
            anigraph.GlDrawDisabledLinks();
        }
        //SampledLoops
        //       if (Nodes0.size()==0)return;

        //ScalarType sizeDiag=anigraph.Mesh().bbox.Diag()*0.0005;
        //        anigraph.GlDrawNodeDebugInfo(Nodes0[0][397],45,20,sizeDiag,true,false,vcg::Color4b::Red);

        //        for (size_t i=0;i<SampledPaths.size();i++)
        //            for (size_t j=0;j<4;j++)
        //                anigraph.GlDrawNodeDebugInfo(SampledPaths[i].nodes[j],45,20,sizeDiag,false,true,vcg::Color4b::Blue);
        if (DrawCandidates)
        {
            for (size_t i=0;i<SampledLoops.size();i++)
            {
                vcg::Color4b Col=vcg::Color4b::ColorRamp(0,SampledLoops.size(),i);
                anigraph.GlDrawPath(SampledLoops[i].LoopP,Col,5,true);
            }
        }

        //        for (size_t i=0;i<anigraph.NumNodes();i++)
        //        {
        //            if (anigraph.IsValid(i))continue;
        //            anigraph.GlDrawNodeDebugInfo(i,45,20,sizeDiag,true,false,vcg::Color4b::Red);
        //        }

        for (size_t i=0;i<ConcaveLoops.size();i++)
            anigraph.GlDrawPath(ConcaveLoops[i],vcg::Color4b::Blue,20,true);

        //        for (size_t i=0;i<GuardPos.size();i++)
        //        {
        //            if (IsGuardConcave[i])
        //                vcg::glColor(vcg::Color4b(0,255,255,255));
        //            else
        //                vcg::glColor(vcg::Color4b(255,0,0,255));

        //            GLDrawFeature(GuardPos[i]);
        //        }
    }
#endif

    void GetSharpNodeSet(std::set<size_t> &SnapNodeSet)
    {
        assert(Nodes0.size()==Nodes1.size());

        for (size_t i=0;i<Nodes0.size();i++)
        {
            for (size_t j=0;j<Nodes0[i].size();j++)
            {
                std::vector<int> TNodes;
                LoopFunctions<MeshType>::GetTangentNodes(anigraph,Nodes0[i][j],TNodes);
                SnapNodeSet.insert(TNodes.begin(),TNodes.end());
            }
        }
        for (size_t i=0;i<Nodes1.size();i++)
            for (size_t j=0;j<Nodes1[i].size();j++)
            {
                std::vector<int> TNodes;
                LoopFunctions<MeshType>::GetTangentNodes(anigraph,Nodes1[i][j],TNodes);
                SnapNodeSet.insert(TNodes.begin(),TNodes.end());
            }
    }

    void GetConcaveEndVert(std::set<size_t> &SharpEndSet)
    {
        SelectConcaveSharpEndPoints();

        //save the set of sherp end vertices
        for (size_t i=0;i<anigraph.Mesh().vert.size();i++)
        {
            if (!anigraph.Mesh().vert[i].IsS())continue;
            SharpEndSet.insert(i);
        }
    }

    void GetConvexEndVert(std::set<size_t> &SharpEndSet)
    {
        size_t NumSel=anigraph.Mesh().SelectEndSharpVert(ConvexEdge);
        (void)NumSel;
        //save the set of sherp end vertices
        for (size_t i=0;i<anigraph.Mesh().vert.size();i++)
        {
            if (!anigraph.Mesh().vert[i].IsS())continue;
            SharpEndSet.insert(i);
        }
    }

    void GetConcaveVert(std::set<size_t> &SharpVertSet)
    {
        SelectConcaveSharpPoints();
        //then the set of vertices on sharp concave features
        for (size_t i=0;i<anigraph.Mesh().vert.size();i++)
        {
            if (!anigraph.Mesh().vert[i].IsS())continue;
            SharpVertSet.insert(i);
        }
    }


    CoordType SnapPos(const size_t IndexNode,
                      std::set<size_t> &SharpVertSet,
                      std::set<size_t> &SharpVertEnd)
    {
        //get the idex of the face
        size_t faceID,M4Dir,edgeIndex;

        anigraph.GetFaceEdgeDir(IndexNode,faceID,edgeIndex,M4Dir);

        int IndexV0=vcg::tri::Index(anigraph.Mesh(),anigraph.Mesh().face[faceID].V0(edgeIndex));
        int IndexV1=vcg::tri::Index(anigraph.Mesh(),anigraph.Mesh().face[faceID].V1(edgeIndex));
        CoordType Pos0=anigraph.Mesh().vert[IndexV0].P();
        CoordType Pos1=anigraph.Mesh().vert[IndexV1].P();

        bool SharpEnd0=(SharpVertEnd.count(IndexV0)>0);
        bool SharpEnd1=(SharpVertEnd.count(IndexV1)>0);
        assert(!((SharpEnd0)&&(SharpEnd1)));

        if (SharpEnd0)return Pos0;
        if (SharpEnd1)return Pos1;

        bool SharpV0=(SharpVertSet.count(IndexV0)>0);
        bool SharpV1=(SharpVertSet.count(IndexV1)>0);
        assert(!(SharpV0 && SharpV1));
        assert((SharpV0 || SharpV1));

        if (SharpV0)return Pos0;
        assert(SharpV1);
        return Pos1;
    }


    void SnapLoops(std::set<size_t> &SnappedNodes)
    {
        //assert(0);
        SnappedNodes.clear();
        //anigraph.Mesh().InitFaceEdgeSelFromFeatureGeo();

        //sharp.InitFaceEdgeSelFromFeatureGeo();

        //        std::set<size_t> ConvexEndSet;
        //        GetConvexEndVert(ConvexEndSet);

        //get concave end vertices
        std::set<size_t> SharpEndSet;
        GetConcaveEndVert(SharpEndSet);

        //then the set of vertices on sharp concave features
        std::set<size_t> SharpVertSet;
        GetConcaveVert(SharpVertSet);

        //then the set of nodes that possibly should be snapped
        std::set<size_t> SnapNodeSet;
        GetSharpNodeSet(SnapNodeSet);

        //then save the ones that must be snapped
        std::vector<std::vector<bool> > MustSnap(ConcaveLoops.size());
        std::vector<std::vector<bool> > HasEndVert(ConcaveLoops.size());
        for (size_t i=0;i<ConcaveLoops.size();i++)
        {
            size_t NumNodes=ConcaveLoops[i].nodes.size();
            MustSnap[i].resize(NumNodes,false);
            HasEndVert[i].resize(NumNodes,false);
            for (size_t j=0;j<NumNodes;j++)
            {
                size_t IndexN=ConcaveLoops[i].nodes[j];

                //get the face and the edge
                size_t IndexF,IndexE,M4Ind;
                anigraph.GetFaceEdgeDir(IndexN,IndexF,IndexE,M4Ind);

                //get the two vertices
                VertexType *v0=anigraph.Mesh().face[IndexF].V0(IndexE);
                VertexType *v1=anigraph.Mesh().face[IndexF].V1(IndexE);
                size_t IndexV0=vcg::tri::Index(anigraph.Mesh(),v0);
                size_t IndexV1=vcg::tri::Index(anigraph.Mesh(),v1);

                //see if close to an end point
                HasEndVert[i][j]=((SharpEndSet.count(IndexV0)>0) || (SharpEndSet.count(IndexV1)>0));

                //check if end pos and within the set
                bool SnapByNode=(SnapNodeSet.count(IndexN)>0);
                if (!SnapByNode)continue;
                MustSnap[i][j]=true;

            }
        }

        //propagate
        bool changed=false;
        do
        {
            changed=false;
            std::cout<<"Step propagate"<<std::endl;
            for (size_t i=0;i<MustSnap.size();i++)
                for (size_t j=0;j<MustSnap[i].size();j++)
                {
                    if (!HasEndVert[i][j])continue;
                    if (MustSnap[i][j])continue;
                    size_t Next=(j+1)%MustSnap[i].size();
                    size_t Prev=(j+MustSnap[i].size()-1)%MustSnap[i].size();
                    if (HasEndVert[i][Next] && MustSnap[i][Next]){MustSnap[i][j]=true;changed=true;}
                    if (HasEndVert[i][Prev] && MustSnap[i][Prev]){MustSnap[i][j]=true;changed=true;}
                }
        }
        while (changed);

        //finally snap
        for (size_t i=0;i<ConcaveLoops.size();i++)
        {
            for (size_t j=0;j<ConcaveLoops[i].nodes.size();j++)
            {
                if (!MustSnap[i][j])continue;
                size_t IndexN=ConcaveLoops[i].nodes[j];
                anigraph.Graph[IndexN]->pos=SnapPos(IndexN,SharpVertSet,SharpEndSet);
                SnappedNodes.insert(IndexN);
            }
        }

        //then snap the convex traced ones
        if ((ConvexCompletionPaths.size()>0)||(ConvexCompletionLoops.size()>0))
        {
            //size_t NumSel=anigraph.Mesh().SelectEndSharpVert();
            size_t NumSel=sharp.SelectEndSharpVert();

            //size_t NumSel=anigraph.Mesh().SelectValence4EndPos();
            assert(NumSel>0);
            for (size_t i=0;i<ConvexCompletionPaths.size();i++)
            {
                for (size_t j=0;j<ConvexCompletionPaths[i].nodes.size();j++)
                {
                    int CurrNode=ConvexCompletionPaths[i].nodes[j];
                    size_t FaceI,EdgeI,M4Dir;
                    anigraph.GetFaceEdgeDir(CurrNode,FaceI,EdgeI,M4Dir);
                    VertexType *v0=anigraph.Mesh().face[FaceI].V0(EdgeI);
                    VertexType *v1=anigraph.Mesh().face[FaceI].V1(EdgeI);
                    if ((!v0->IsS()) && (!v1->IsS()))continue;
                    assert(!(v0->IsS() && v1->IsS()));
                    if (v0->IsS())anigraph.Graph[CurrNode]->pos=v0->P();
                    if (v1->IsS())anigraph.Graph[CurrNode]->pos=v1->P();
                }
            }
            for (size_t i=0;i<ConvexCompletionLoops.size();i++)
            {
                for (size_t j=0;j<ConvexCompletionLoops[i].nodes.size();j++)
                {
                    int CurrNode=ConvexCompletionLoops[i].nodes[j];
                    size_t FaceI,EdgeI,M4Dir;
                    anigraph.GetFaceEdgeDir(CurrNode,FaceI,EdgeI,M4Dir);
                    VertexType *v0=anigraph.Mesh().face[FaceI].V0(EdgeI);
                    VertexType *v1=anigraph.Mesh().face[FaceI].V1(EdgeI);
                    if ((!v0->IsS()) && (!v1->IsS()))continue;
                    assert(!(v0->IsS() && v1->IsS()));
                    if (v0->IsS()){anigraph.Graph[CurrNode]->pos=v0->P();SnappedNodes.insert(CurrNode);}
                    if (v1->IsS()){anigraph.Graph[CurrNode]->pos=v1->P();SnappedNodes.insert(CurrNode);}
                }
            }
        }

    }

//    std::vector<CoordType> PD1C,PD2C;
//    std::vector<ScalarType> K1,K2;

//    void InitCurvVectors()
//    {
//        //typedef typename MeshType::ScalarType ScalarType;
//        typedef typename MeshType::CoordType CoordType;

//        //first compute the curvature on the mesh
//        Eigen::MatrixXi F;
//        Eigen::MatrixXd V;
//        Eigen::MatrixXd PD1,PD2,PV1,PV2;
//        MeshToMatrix<MeshType>::GetTriMeshData(anigraph.Mesh(),F,V);
//        igl::principal_curvature(V,F,PD1,PD2,PV1,PV2,3);

//        PD1C.clear();
//        PD2C.clear();
//        K1.clear();
//        K2.clear();
//        for (int i=0;i<PD1.rows();i++)
//        {
//            PD1C.push_back(CoordType(PD1(i,0),PD1(i,1),PD1(i,2)));
//            PD2C.push_back(CoordType(PD2(i,0),PD2(i,1),PD2(i,2)));
//            K1.push_back(PV1(i,0));
//            K2.push_back(PV2(i,0));
//        }
//    }

    void GetLoopsGreedyAdd(bool _AvoidCrossIntersections,
                           bool _CloseConvexEndPoints,
                           bool _CloseOnlyWithOrtho,
                           std::vector<PathType> &_ConcaveLoops,
                           std::vector<PathType> &_ProblematicLoops,
                           std::vector<PathType> &_ConvexPath,
                           std::vector<std::vector<vcg::face::Pos<FaceType> > > &_ClosedConcavePos,
                           std::vector<std::vector<vcg::face::Pos<FaceType> > > &_UnsampledConcavePos,
                           std::vector<std::vector<vcg::face::Pos<FaceType> > > &_ConvexPos,
                           ScalarType _AlignFactor,
                           ScalarType _MaxAngle)
    {
//        InitCurvVectors();

        MaxAngle=_MaxAngle;
        AlignFactor=_AlignFactor;
        AvoidCrossIntersections=_AvoidCrossIntersections;
        CloseConvexEndPoints=_CloseConvexEndPoints;
        CloseOnlyWithOrtho=_CloseOnlyWithOrtho;

        _ConcaveLoops.clear();
        _ClosedConcavePos.clear();
        _UnsampledConcavePos.clear();
        _ConvexPos.clear();
        _ProblematicLoops.clear();


        std::cout<<"Initializing Graph"<<std::endl;

        InitStructures();


        TimeTrace=0;
        TimeCheck=0;


        TraceLoops();

        //ColorSharpByCovering();

        _ProblematicLoops=ProblematicLoops;
        _ConcaveLoops=ConcaveLoops;

        //then set the unsampled
        for (size_t i=0;i<Covered0.size();i++)
        {
            if ((Covered0[i]) || (Covered1[i]))continue;
            _UnsampledConcavePos.push_back(ConcavePos[i]);
        }
        for (size_t i=0;i<GuardPos.size();i++)
        {
            if (IsGuardConcave[i])
                _ClosedConcavePos.push_back(GuardPos[i]);
            else
                _ConvexPos.push_back(GuardPos[i]);
        }
        if (CloseConvexEndPoints)
        {
            ConvexClose();
            _ConcaveLoops.insert(_ConcaveLoops.end(),ConvexCompletionLoops.begin(),ConvexCompletionLoops.end());
            _ConvexPath.insert(_ConvexPath.end(),ConvexCompletionPaths.begin(),ConvexCompletionPaths.end());
        }

        std::cout<<"*** FINAL STATS Sharp Tracing***"<<std::endl;
        std::cout<<"* Closed Loops From Begin "<<_ClosedConcavePos.size()<<std::endl;
        std::cout<<"* Loops Traced "<<_ConcaveLoops.size()<<std::endl;
        std::cout<<"* Concave UnTraced "<<_UnsampledConcavePos.size()<<std::endl;
        std::cout<<"* Convex Features "<<_ConvexPos.size()<<std::endl;
    }

};

#endif //LOOP_OPTIMIZER
