#ifndef LOOP_COMMON_FUNCTIONS
#define LOOP_COMMON_FUNCTIONS

#include <tracing_field/graph/anisotropic_graph.h>
//#include <tracing_field/anisotropic_geodesic.h>
#include <tracing_field/conflict_finder.h>
#include <wrap/io_trimesh/export.h>
#include <unordered_map>

enum SharpFeatureType{ConvexEdge,ConcaveEdge,NoSharp};

/**
 * @brief The class that performs parametrization by considering the connections within the graph
 */
template < class MeshType >
class LoopFunctions {

public:

    typedef typename MeshType::CoordType    CoordType;
    typedef typename MeshType::VertexType   VertexType;
    typedef typename MeshType::FaceType     FaceType;
    typedef typename MeshType::ScalarType   ScalarType;
    typedef typename AnisotropicGraph<MeshType>::Path PathType;
    typedef typename AnisotropicGraph<MeshType>::NeighborsIterator NeighborsIterator;
    typedef typename AnisotropicGraph<MeshType>::NodeListIterator NodeListIterator;
    typedef typename AnisotropicGraph<MeshType>::NodeListConstIterator NodeListConstIterator;

    typedef typename ConflictFinder<MeshType>::ConflictInfo ConflictInfo;

    template <class EdgeMeshType>
    static void SetEdgeSelFromEMesh(MeshType &tri_mesh,
                                    EdgeMeshType &EMesh)
    {
        vcg::tri::UpdateFlags<MeshType>::FaceClearFaceEdgeS(tri_mesh);
        std::set<std::pair<CoordType,CoordType> > EdgeSet;
        for (size_t i=0;i<EMesh.edge.size();i++)
        {
            CoordType P0=EMesh.edge[i].P(0);
            CoordType P1=EMesh.edge[i].P(1);
            EdgeSet.insert(std::pair<CoordType,CoordType>(std::min(P0,P1),std::max(P0,P1)));
        }
        for (size_t i=0;i<tri_mesh.face.size();i++)
            for (size_t j=0;j<3;j++)
            {
                CoordType P0=EMesh.edge[i].P(0);
                CoordType P1=EMesh.edge[i].P(1);
                std::pair<CoordType,CoordType> key(std::min(P0,P1),std::max(P0,P1));
                if (EdgeSet.count(key)==0)continue;
                tri_mesh.face[i].SetFaceEdgeS(j);
            }
    }

    static void RetrieveGuardNodes(const AnisotropicGraph<MeshType> &anigraph,
                                   int currF,int currE,std::vector<size_t> &nodesI)
    {
        CoordType Pos0=anigraph.Mesh().face[currF].P0(currE);
        CoordType Pos1=anigraph.Mesh().face[currF].P1(currE);
        CoordType currDir=Pos0-Pos1;
        currDir.Normalize();
        std::vector<size_t> nodesI0,nodesI1;
        anigraph.faceEdgeNodesDirCoord(currF,currE,currDir,nodesI0);
        anigraph.faceEdgeNodesDirCoord(currF,currE,-currDir,nodesI1);
        nodesI=std::vector<size_t>(nodesI0.begin(),nodesI0.end());
        nodesI.insert(nodesI.end(),nodesI1.begin(),nodesI1.end());
    }

    static void RetrieveCrossNodes(const AnisotropicGraph<MeshType> &anigraph,
                                   int currF,int currE,std::vector<size_t> &nodesI)
    {
        CoordType Pos0=anigraph.Mesh().face[currF].P0(currE);
        CoordType Pos1=anigraph.Mesh().face[currF].P1(currE);
        CoordType currDir=Pos0-Pos1;
        currDir.Normalize();

        //then rotate
        CoordType Norm=anigraph.Mesh().face[currF].N();
        currDir=currDir^Norm;
        currDir.Normalize();

        std::vector<size_t> nodesI0,nodesI1;
        anigraph.faceEdgeNodesDirCoord(currF,currE,currDir,nodesI0);
        anigraph.faceEdgeNodesDirCoord(currF,currE,-currDir,nodesI1);
        nodesI=std::vector<size_t>(nodesI0.begin(),nodesI0.end());
        nodesI.insert(nodesI.end(),nodesI1.begin(),nodesI1.end());
    }

    static void RetrieveGuardNodes(const AnisotropicGraph<MeshType> &anigraph,
                                   const std::vector<vcg::face::Pos<FaceType> > &FeaturePos,
                                   std::vector<size_t> &GuardNodes)
    {

        for (size_t i=0;i<FeaturePos.size();i++)
        {
            vcg::face::Pos<FaceType> CurrP=FeaturePos[i];
            int currF=vcg::tri::Index(anigraph.Mesh(),CurrP.F());
            int currE=CurrP.E();
            std::vector<size_t> nodesI0,nodesI1;
            RetrieveGuardNodes(anigraph,currF,currE,nodesI0);
            //then go on the opposite side
            int OppFace=vcg::tri::Index(anigraph.Mesh(),anigraph.Mesh().face[currF].FFp(currE));
            int OppE=anigraph.Mesh().face[currF].FFi(currE);
            RetrieveGuardNodes(anigraph,OppFace,OppE,nodesI1);
            GuardNodes.insert(GuardNodes.end(),nodesI0.begin(),nodesI0.end());
            GuardNodes.insert(GuardNodes.end(),nodesI1.begin(),nodesI1.end());
        }
    }

    static void RetrieveCrossNodes(const AnisotropicGraph<MeshType> &anigraph,
                                   const std::vector<vcg::face::Pos<FaceType> > &FeaturePos,
                                   std::vector<size_t> &GuardNodes)
    {

        for (size_t i=0;i<FeaturePos.size();i++)
        {
            vcg::face::Pos<FaceType> CurrP=FeaturePos[i];
            int currF=vcg::tri::Index(anigraph.Mesh(),CurrP.F());
            int currE=CurrP.E();
            std::vector<size_t> nodesI0,nodesI1;
            RetrieveCrossNodes(anigraph,currF,currE,nodesI0);
            //then go on the opposite side
            int OppFace=vcg::tri::Index(anigraph.Mesh(),anigraph.Mesh().face[currF].FFp(currE));
            int OppE=anigraph.Mesh().face[currF].FFi(currE);
            RetrieveCrossNodes(anigraph,OppFace,OppE,nodesI1);
            GuardNodes.insert(GuardNodes.end(),nodesI0.begin(),nodesI0.end());
            GuardNodes.insert(GuardNodes.end(),nodesI1.begin(),nodesI1.end());
        }
    }

    //    static void SetBarrierSharp(AnisotropicGraph<MeshType> &anigraph,
    //                                const std::vector<std::vector<vcg::face::Pos<FaceType> > > &FeaturePos,
    //                                const std::vector<size_t> &ExcludeNodes)
    //    {
    //        //collect all guard nodes
    //        //std::vector<size_t> GuardNodes;
    //        std::set<std::pair<size_t,size_t> > FaceDir;
    //        for (size_t i=0;i<FeaturePos.size();i++)
    //            for (size_t j=0;j<FeaturePos[i].size();j++)
    //        {
    //            size_t FaceId=vcg::tri::Index(anigraph.Mesh(),FeaturePos[i][j].F());
    //            int M4Dir=getFaceEdgeOrientedDir(FeaturePos[i][j]);
    //            FaceDir.insert(std::pair<size_t,size_t>(FaceId,M4Dir));
    //            FaceDir.insert(std::pair<size_t,size_t>(FaceId,(M4Dir+2)%4));
    //        }
    //            //RetrieveGuardNodes(anigraph,FeaturePos[i],GuardNodes);

    ////        std::set<std::pair<size_t,size_t> > FaceDir;
    ////        for (size_t i=0;i<GuardNodes.size();i++)
    ////        {
    ////            int IndexNode=GuardNodes[i];
    ////            assert(anigraph.IsEdgeNode(IndexNode));
    ////            size_t FaceId;
    ////            size_t M4Dir;
    ////            anigraph.GetFaceDir(IndexNode,FaceId,M4Dir);
    ////            FaceDir.insert(std::pair<size_t,size_t>(FaceId,M4Dir));
    ////            FaceDir.insert(std::pair<size_t,size_t>(FaceId,(M4Dir+2)%4));
    ////        }

    //        BarrierPaths(anigraph,FaceDir,ExcludeNodes);
    //    }

        static void DisableFaceDir(AnisotropicGraph<MeshType> &anigraph,
                                 const std::set<std::pair<size_t,size_t> > &FaceDir)

        {
            //then invalidate all nodes considering such directions
            for (size_t i=0;i<anigraph.NumNodes();i++)
            {
                size_t FaceId;
                size_t M4Dir;
                if (!anigraph.IsValid(i))continue;
                if (!anigraph.IsEdgeNode(i))continue;
                anigraph.GetFaceDir(i,FaceId,M4Dir);
                std::pair<size_t,size_t> key(FaceId,M4Dir);
                if (FaceDir.count(key)==0)continue;
                std::vector<int> TangNodes;
                GetTangentNodes(anigraph,i,TangNodes);
                for (size_t j=0;j<TangNodes.size();j++)
                    anigraph.SetValid(TangNodes[j],false);
//                if (To_Exclude.count(i)==0)
//                {
//                    if (Nodes!=NULL)Nodes->push_back(i);

//                    anigraph.SetValid(i,false);
//                }
//                int IndexOpp=anigraph.OppositeNode(i);
//                assert(anigraph.IsEdgeNode(IndexOpp));

//                if (To_Exclude.count(IndexOpp)==0)
//                {
//                    if (Nodes!=NULL)Nodes->push_back(IndexOpp);
//                    anigraph.SetValid(IndexOpp,false);
//                }
            }
            anigraph.InvalidateArcsOnNonValidNodes();
        }

    static void DisableFaceSharp(AnisotropicGraph<MeshType> &anigraph,
                                const std::vector<std::vector<vcg::face::Pos<FaceType> > > &FeaturePos)
    {
        //collect all guard nodes
        //std::vector<size_t> GuardNodes;
        std::set<std::pair<size_t,size_t> > FaceDir;
        for (size_t i=0;i<FeaturePos.size();i++)
            for (size_t j=0;j<FeaturePos[i].size();j++)
        {
            size_t FaceId=vcg::tri::Index(anigraph.Mesh(),FeaturePos[i][j].F());
            int M4Dir=getFaceEdgeOrientedDir(FeaturePos[i][j]);
            FaceDir.insert(std::pair<size_t,size_t>(FaceId,M4Dir));
            FaceDir.insert(std::pair<size_t,size_t>(FaceId,(M4Dir+2)%4));
        }

        DisableFaceDir(anigraph,FaceDir);
    }

    //    struct FaceEdgeDir
    //    {
    //      size_t IndexF;
    //      size_t IndexE;
    //      size_t M4Dir;
    //      inline bool operator < (const FaceEdgeDir &FEdir)const
    //      {
    //          if ((IndexF==FEdir.IndexF)&&(IndexE==FEdir.IndexE))
    //              return (M4Dir<FEdir.M4Dir);
    //          if ((IndexF==FEdir.IndexF))
    //              return (IndexE<FEdir.IndexE);
    //          return (IndexF<FEdir.IndexF);
    //      }

    //      inline bool operator == (const FaceEdgeDir &FEdir)const
    //      {
    //          return ((IndexF==FEdir.IndexF)&&
    //                  (IndexE==FEdir.IndexE)&&
    //                  (M4Dir==FEdir.M4Dir));
    //      }

    //      FaceEdgeDir(size_t &_IndexF,size_t &_IndexE,size_t &_M4Dir)
    //      {
    //          IndexF=_IndexF;
    //          IndexE=_IndexE;
    //          M4Dir=_M4Dir;
    //      }
    //    };

    static void SetBarrierSharp(AnisotropicGraph<MeshType> &anigraph,
                                const std::vector<std::vector<vcg::face::Pos<FaceType> > > &FeaturePos,
                                const std::vector<size_t> &ExcludeNodes)
    {
        std::set<size_t> ExcludeSet(ExcludeNodes.begin(),ExcludeNodes.end());
        std::vector<size_t> BoundaryNodes;
        for (size_t i=0;i<FeaturePos.size();i++)
            for (size_t j=0;j<FeaturePos[i].size();j++)
            {
                int FaceId=vcg::tri::Index(anigraph.Mesh(),FeaturePos[i][j].F());
                int M4Dir=getFaceEdgeOrientedDir(FeaturePos[i][j]);
                int IndexE=FeaturePos[i][j].E();
                anigraph.faceEdgeNodesDir(FaceId,IndexE,M4Dir,BoundaryNodes);
                M4Dir=(M4Dir+2)%4;
                anigraph.faceEdgeNodesDir(FaceId,IndexE,M4Dir,BoundaryNodes);

                //go to the other side
                int FaceId1=vcg::tri::Index(anigraph.Mesh(),anigraph.Mesh().face[FaceId].FFp(IndexE));
                int IndexE1=anigraph.Mesh().face[FaceId].FFi(IndexE);
                if (IndexE1<0)continue;
                vcg::face::Pos<FaceType> PosOpp(&anigraph.Mesh().face[FaceId1],IndexE1);
                size_t M4Dir1=getFaceEdgeOrientedDir(PosOpp);
                anigraph.faceEdgeNodesDir(FaceId1,IndexE1,M4Dir1,BoundaryNodes);
                M4Dir1=(M4Dir1+2)%4;
                anigraph.faceEdgeNodesDir(FaceId1,IndexE1,M4Dir1,BoundaryNodes);
                //            size_t IndexE=FeaturePos[i].E();
                //            ExcludeDomain.push_back(FaceEdgeDir(FaceId,IndexE,M4Dir));
            }
        std::sort(BoundaryNodes.begin(),BoundaryNodes.end());
        std::vector<size_t>::iterator iteTerm=std::unique(BoundaryNodes.begin(),BoundaryNodes.end());
        BoundaryNodes.erase(iteTerm, BoundaryNodes.end());

        for (size_t i=0;i<BoundaryNodes.size();i++)
        {
            int IndexNode=BoundaryNodes[i];
            if (ExcludeSet.count(IndexNode)>0)continue;

            anigraph.SetValid(IndexNode,false);
            int OppNode=anigraph.OppositeNode(IndexNode);
            anigraph.SetValid(OppNode,false);
            //            int IndexOpp=anigraph.OppositeNode(i);
            //            assert(anigraph.IsEdgeNode(IndexOpp));

            //            if (To_Exclude.count(IndexOpp)==0)
            //                anigraph.SetValid(IndexOpp,false);
        }
        anigraph.InvalidateArcsOnNonValidNodes();

        //        std::set<std::pair<size_t,size_t> > FaceDir;
        //        for (size_t i=0;i<GuardNodes.size();i++)
        //        {
        //            int IndexNode=GuardNodes[i];
        //            assert(anigraph.IsEdgeNode(IndexNode));
        //            size_t FaceId;
        //            size_t M4Dir;
        //            anigraph.GetFaceDir(IndexNode,FaceId,M4Dir);
        //            FaceDir.insert(std::pair<size_t,size_t>(FaceId,M4Dir));
        //            FaceDir.insert(std::pair<size_t,size_t>(FaceId,(M4Dir+2)%4));
        //        }

        //        BarrierPaths(anigraph,FaceDir,ExcludeNodes);
    }


    static void GetNodesOfFaceDir(AnisotropicGraph<MeshType> &anigraph,
                                  const std::set<std::pair<size_t,size_t> > &FaceDir,
                                  std::vector<size_t> &Nodes)
    {
        for (size_t i=0;i<anigraph.NumNodes();i++)
        {
            size_t FaceId;
            size_t M4Dir;
            //            if (!anigraph.IsValid(i))continue;
            if (!anigraph.IsEdgeNode(i))continue;
            anigraph.GetFaceDir(i,FaceId,M4Dir);
            std::pair<size_t,size_t> key(FaceId,M4Dir);
            if (FaceDir.count(key)==0)continue;
            int IndexOpp=anigraph.OppositeNode(i);
            assert(anigraph.IsEdgeNode(IndexOpp));
            Nodes.push_back(i);
            Nodes.push_back(IndexOpp);
        }
    }


//    static void BarrierPaths(AnisotropicGraph<MeshType> &anigraph,
//                             const std::set<std::pair<size_t,size_t> > &FaceDir,
//                             const std::vector<size_t> &ExcludeNodes,
//                             std::vector<size_t> *Nodes=NULL)

//    {
//        std::set<size_t> To_Exclude(ExcludeNodes.begin(),ExcludeNodes.end());
//        //then invalidate all nodes considering such directions
//        for (size_t i=0;i<anigraph.NumNodes();i++)
//        {
//            size_t FaceId;
//            size_t M4Dir;
//            if (!anigraph.IsValid(i))continue;
//            if (!anigraph.IsEdgeNode(i))continue;
//            anigraph.GetFaceDir(i,FaceId,M4Dir);
//            std::pair<size_t,size_t> key(FaceId,M4Dir);
//            if (FaceDir.count(key)==0)continue;

//            if (To_Exclude.count(i)==0)
//            {
//                if (Nodes!=NULL)Nodes->push_back(i);

//                anigraph.SetValid(i,false);
//            }
//            int IndexOpp=anigraph.OppositeNode(i);
//            assert(anigraph.IsEdgeNode(IndexOpp));

//            if (To_Exclude.count(IndexOpp)==0)
//            {
//                if (Nodes!=NULL)Nodes->push_back(IndexOpp);
//                anigraph.SetValid(IndexOpp,false);
//            }
//        }
//        anigraph.InvalidateArcsOnNonValidNodes();
//    }

//    static void BarrierPaths(AnisotropicGraph<MeshType> &anigraph,
//                             const std::vector<PathType> &TestLoops,
//                             const std::vector<bool> &IsLoop,
//                             std::vector<size_t> *Nodes=NULL)
//    {
//        std::set<std::pair<size_t,size_t> > FaceDir;
//        for (size_t i=0;i<TestLoops.size();i++)
//        {
//            std::vector<std::pair<size_t,size_t> > CurrFaceDir;
//            anigraph.GetFacesDirM2(TestLoops[i],CurrFaceDir,IsLoop[i]);

//            for (size_t j=0;j<CurrFaceDir.size();j++)
//                FaceDir.insert(CurrFaceDir[j]);
//        }
//        std::vector<size_t> ExclNodes;
//        BarrierPaths(anigraph,FaceDir,ExclNodes,Nodes);
//    }

    static void GetTangentNodes(AnisotropicGraph<MeshType> &anigraph,
                                int IndexNode,
                                std::vector<int> &TangNodes)
    {
        TangNodes.clear();
        assert(IndexNode>=0);
        assert(IndexNode<(int)anigraph.NumNodes());
        int TangTest[4];
        TangTest[0]=IndexNode;
        //std::cout<<"A"<<TangTest[0]<<std::endl;
        TangTest[1]=anigraph.OppositeNode(TangTest[0]);
        //std::cout<<"B"<<TangTest[1]<<std::endl;
        TangTest[2]=anigraph.OppositeNodeSameF(TangTest[0]);
        //std::cout<<"C"<<TangTest[2]<<std::endl;
        TangTest[3]=anigraph.OppositeNodeSameF(TangTest[1]);
        //std::cout<<"D"<<TangTest[3]<<std::endl;

        for (size_t k=0;k<4;k++)
        {
            if (TangTest[k]<0)continue;
            TangNodes.push_back(TangTest[k]);
        }
    }

    static void GetCrossNodes(AnisotropicGraph<MeshType> &anigraph,
                              int IndexNode,
                              std::vector<int> &CrossNodes)
    {
        std::vector<int> TangNodes;
        GetTangentNodes(anigraph,IndexNode,TangNodes);
        std::set<int> TangNodesSet=std::set<int>(TangNodes.begin(),TangNodes.end());

        CrossNodes.clear();
        CoordType pos=anigraph.NodePos(IndexNode);

        NodeListConstIterator First,Last,Curr;
        anigraph.GetCoordNodes(pos,First,Last);
        for (Curr=First;Curr!=Last;Curr++)
        {
            size_t IndexN=(*Curr);
            assert(IndexN>=0);
            assert(IndexN<anigraph.NumNodes());
            if (TangNodesSet.count((int)IndexN)>0)continue;
            CrossNodes.push_back(IndexN);
        }
    }

    static void GetTangentNodes(AnisotropicGraph<MeshType> &anigraph,
                                const PathType &TestLoop,
                                std::vector< std::vector<int> > &TangentNodes)
    {
        TangentNodes.clear();
        for (size_t j=0;j<TestLoop.nodes.size();j++)
        {
            std::vector<int> CurrTangNodes;
            GetTangentNodes(anigraph,TestLoop.nodes[j],CurrTangNodes);
            TangentNodes.push_back(CurrTangNodes);
        }
    }

//    static void SetBarrierForPaths(AnisotropicGraph<MeshType> &anigraph,
//                                   const std::vector<PathType> &BarrierPaths,
//                                   const std::vector<bool> &IsLoop)
//    {
//        //if (Nodes!=NULL)Nodes->clear();

//        //first disable all tangent nodes
//        for (size_t i=0;i<BarrierPaths.size();i++)
//        {
//            std::vector< std::vector<int> > TangentNodes;
//            GetTangentNodes(anigraph,BarrierPaths[i],TangentNodes);
//            for (size_t j=0;j<TangentNodes.size();j++)
//                for (size_t k=0;k<TangentNodes[j].size();k++)
//                    anigraph.SetValid(TangentNodes[j][k],false);
//        }

//        anigraph.InvalidateArcsOnNonValidNodes();
//        //        std::set<std::pair<size_t,size_t> > FaceDir;
//        //        for (size_t i=0;i<TestLoops.size();i++)
//        //        {
//        //            std::vector<std::pair<size_t,size_t> > CurrFaceDir;
//        //            anigraph.GetFacesDirM2(TestLoops[i],CurrFaceDir,IsLoop[i]);

//        //            for (size_t j=0;j<CurrFaceDir.size();j++)
//        //                FaceDir.insert(CurrFaceDir[j]);
//        //        }
//        //        std::vector<size_t> ExclNodes;
//        //        BarrierPaths(anigraph,FaceDir,ExclNodes,Nodes);
//    }

    //    static void SetBarrierForPaths(AnisotropicGraph<MeshType> &anigraph,
    //                                   const std::vector<PathType> &TestLoops,
    //                                   const std::vector<bool> &IsLoop,
    //                                   std::vector<size_t> *Nodes=NULL)
    //    {
    //        if (Nodes!=NULL)Nodes->clear();

    //        //first disable all the nodes
    //        for (size_t i=0;i<TestLoops.size();i++)
    //            for (size_t j=0;j<TestLoops[i].nodes.size();j++)
    //            {
    //                size_t DisableNodes[4];
    //                DisableNodes[0]=TestLoops[i].nodes[j];
    //                DisableNodes[1]=anigraph.OppositeNode(DisableNodes[0]);
    //                DisableNodes[2]=anigraph.OppositeNodeSameF(DisableNodes[0]);
    //                DisableNodes[3]=anigraph.OppositeNodeSameF(DisableNodes[1]);
    //                for (size_t k=0;k<4;k++)
    //                {
    //                    if (DisableNodes[k]<0)continue;
    //                    anigraph.SetValid(DisableNodes[k],false);
    //                }
    //            }

    //        anigraph.InvalidateArcsOnNonValidNodes();
    //        //        std::set<std::pair<size_t,size_t> > FaceDir;
    //        //        for (size_t i=0;i<TestLoops.size();i++)
    //        //        {
    //        //            std::vector<std::pair<size_t,size_t> > CurrFaceDir;
    //        //            anigraph.GetFacesDirM2(TestLoops[i],CurrFaceDir,IsLoop[i]);

    //        //            for (size_t j=0;j<CurrFaceDir.size();j++)
    //        //                FaceDir.insert(CurrFaceDir[j]);
    //        //        }
    //        //        std::vector<size_t> ExclNodes;
    //        //        BarrierPaths(anigraph,FaceDir,ExclNodes,Nodes);
    //    }


    static bool IsClosedLoop(std::vector<vcg::face::Pos<FaceType> > &PosSeq)
    {
        for (size_t i=0;i<PosSeq.size();i++)
        {
            VertexType *v0=PosSeq[i].V();
            VertexType *v1=PosSeq[i].VFlip();
            v0->Q()=0;
            v1->Q()=0;
        }

        for (size_t i=0;i<PosSeq.size();i++)
        {
            VertexType *v0=PosSeq[i].V();
            VertexType *v1=PosSeq[i].VFlip();
            v0->Q()+=1;
            v1->Q()+=1;
        }

        for (size_t i=0;i<PosSeq.size();i++)
        {
            VertexType *v0=PosSeq[i].V();
            VertexType *v1=PosSeq[i].VFlip();
            if (v0->Q()==1)return false;
            if (v1->Q()==1)return false;
        }
        return true;
    }

    static bool IsClosedLoop(MeshType &m,
                             const std::vector<std::pair<size_t,size_t> > &EdgeFaceSeq)
    {
        std::vector<vcg::face::Pos<FaceType> > PosSeq;
        for (size_t i=0;i<EdgeFaceSeq.size();i++)
        {
            size_t IndexF=EdgeFaceSeq[i].first;
            size_t IndexE=EdgeFaceSeq[i].second;
            PosSeq.push_back(vcg::face::Pos<FaceType>(&m.face[IndexF],IndexE));
        }
        return IsClosedLoop(PosSeq);
    }
    //    static void SortSharpFeature(MeshType Mesh,std::vector<vcg::face::Pos<FaceType> > &PosSeq)
    //    {
    //       if (!IsClosedLoop(PosSeq))
    //       {

    //           std::vector<size_t> SharpEndP;
    //           GetSharpEndPoints(Mesh,PosSeq,SharpEndP);
    //       }
    //    }

    static ScalarType PosLenght(std::vector<vcg::face::Pos<FaceType> > &PosSeq)
    {
        ScalarType CurrL=0;
        for (size_t i=0;i<PosSeq.size();i++)
        {
            VertexType *v0=PosSeq[i].V();
            VertexType *v1=PosSeq[i].VFlip();
            CurrL+=(v0->P()-v1->P()).Norm();
        }
        return CurrL;
    }

    static void FollowCrease(MeshType &Mesh,
                             const vcg::face::Pos<FaceType> &StartPos,
                             std::vector<size_t> &FacesI)
    {
        FacesI.clear();

        //the first one should be crease
        vcg::face::Pos<FaceType> Pos=StartPos;
        //assert(Pos.F()->IsCrease(Pos.E()));
        assert(Pos.F()->IsFaceEdgeS(Pos.E()));

        //set the first index
        bool has_terminated=false;
        do{
            size_t FIndex=vcg::tri::Index(Mesh,Pos.F());
            FacesI.push_back(FIndex);

            FaceType *f0=Pos.F();

            //go to next one
            //Pos.NextCrease();
            Pos.NextEdgeS();

            //in case of a cycle terminated when came back
            has_terminated|=(Pos==StartPos);

            //or it found the end of the feature
            FaceType *f1=Pos.FFlip();
            has_terminated|=(f0==f1);
        }while(!has_terminated);
    }

    static void GetFaceSharpPartitionSide(MeshType &Mesh,
                                          const vcg::face::Pos<FaceType> &SharpPos,
                                          std::vector<size_t> &IndexF)
    {
        vcg::face::Pos<FaceType> Pos0=SharpPos;
        vcg::face::Pos<FaceType> Pos1=Pos0;
        Pos1.FlipV();

        std::vector<size_t> IndexF0,IndexF1;

        //go on one side
        FollowCrease(Mesh,Pos0,IndexF0);

        //and on the other side
        FollowCrease(Mesh,Pos1,IndexF1);

        IndexF=IndexF0;
        IndexF.insert(IndexF.end(),IndexF1.begin(),IndexF1.end());
    }

    //set the sharp edge as creases
    static void InitSharpFeaturesAsCreases(MeshType &Mesh,
                                           const std::vector<vcg::face::Pos<FaceType> > &FeaturePos)
    {
        //vcg::tri::UpdateFlags<MeshType>::FaceClearCreases(Mesh);
        vcg::tri::UpdateFlags<MeshType>::FaceClearFaceEdgeS(Mesh);

        for (size_t j=0;j<FeaturePos.size();j++)
        {

            vcg::face::Pos<FaceType> currPos=FeaturePos[j];
            size_t F0=vcg::tri::Index(Mesh,currPos.F());
            size_t E0=currPos.E();

            //go on the other side
            currPos.FlipF();
            size_t F1=vcg::tri::Index(Mesh,currPos.F());
            size_t E1=currPos.E();

            //set faux edges
            //Mesh.face[F0].SetCrease(E0);
            //Mesh.face[F1].SetCrease(E1);
            Mesh.face[F0].SetFaceEdgeS(E0);
            Mesh.face[F1].SetFaceEdgeS(E1);
        }
    }

    static void GetSharpEndPoints(MeshType &Mesh,
                                  const std::vector<vcg::face::Pos<FaceType> > &SharpPos,
                                  std::vector<size_t> &SharpEndP)
    {
        SharpEndP.clear();
        for (size_t i=0;i<SharpPos.size();i++)
        {
            VertexType *v0=SharpPos[i].V();
            VertexType *v1=SharpPos[i].VFlip();
            v0->Q()=0;
            v1->Q()=0;
        }

        for (size_t i=0;i<SharpPos.size();i++)
        {
            VertexType *v0=SharpPos[i].V();
            VertexType *v1=SharpPos[i].VFlip();
            v0->Q()+=1;
            v1->Q()+=1;
        }

        for (size_t i=0;i<SharpPos.size();i++)
        {
            VertexType *v0=SharpPos[i].V();
            VertexType *v1=SharpPos[i].VFlip();
            if (v0->Q()==1)SharpEndP.push_back(vcg::tri::Index(Mesh,v0));
            if (v1->Q()==1)SharpEndP.push_back(vcg::tri::Index(Mesh,v1));
        }
    }

    static void SelectSharpVertices(MeshType &Mesh,
                                    const std::vector<vcg::face::Pos<FaceType> > &SharpPos)
    {
        vcg::tri::UpdateSelection<MeshType>::VertexClear(Mesh);
        for (size_t i=0;i<SharpPos.size();i++)
        {
            VertexType *v0=SharpPos[i].V();
            VertexType *v1=SharpPos[i].VFlip();
            v0->SetS();
            v1->SetS();
        }
    }

    static void PropagatePropertiesOnSharpPartition(MeshType &Mesh,
                                                    const std::vector<vcg::face::Pos<FaceType> > &SharpPos,
                                                    const std::vector<size_t> SharpEndP,
                                                    std::vector<int> &PreferredDir,
                                                    std::vector<size_t> &IndexF)
    {

        //select sharp feature vertices
        SelectSharpVertices(Mesh,SharpPos);

        //        //get end points
        //        std::vector<size_t> SharpEndP;
        //        GetSharpEndPoints(Mesh,SharpPos,SharpEndP);

        //remove end points from selection
        for(size_t i=0;i<SharpEndP.size();i++)
            Mesh.vert[SharpEndP[i]].ClearS();

        //set one to quality to the initialized ones
        vcg::tri::UpdateSelection<MeshType>::FaceClear(Mesh);
        for (size_t i=0;i<IndexF.size();i++)
            Mesh.face[IndexF[i]].SetS();

        //set as crease the sharp features
        InitSharpFeaturesAsCreases(Mesh,SharpPos);

        bool HasMoved=false;
        do
        {
            HasMoved=false;

            for (size_t j=0;j<Mesh.face.size();j++)
            {
                FaceType *CurrF=&Mesh.face[j];

                //already tagged
                if (CurrF->IsS())continue;

                //get the neighbours
                for (size_t k=0;k<3;k++)
                {
                    //if crease then cannot pass on other side
                    //if (CurrF->IsCrease(k))continue;
                    if (CurrF->IsFaceEdgeS(k))continue;

                    FaceType *NextF=CurrF->FFp(k);
                    if (!NextF->IsS())continue;

                    //see if has one vertex on the sharp feature
                    bool HasV=false;

                    VertexType *v0=Mesh.face[j].V0(k);
                    VertexType *v1=Mesh.face[j].V1(k);

                    //all but not at the end of sharp features
                    HasV|=v0->IsS();
                    HasV|=v1->IsS();


                    //if it has such vertex than change the partition
                    if (!HasV)continue;

                    //otherwise propagate

                    //get preferred dir
                    int NextIndex=vcg::tri::Index(Mesh,NextF);



                    int M4DirNextI=PreferredDir[NextIndex];
                    assert(M4DirNextI>=0);
                    assert(M4DirNextI<=4);

                    CoordType M4DirNext=vcg::tri::CrossField<MeshType>::CrossVector(*NextF,M4DirNextI);

                    //transport preferred dir to next face
                    int M4DirCurr=vcg::tri::CrossField<MeshType>::FollowDirectionI(*NextF,*CurrF,M4DirNext);
                    assert(M4DirCurr>=0);
                    assert(M4DirCurr<4);
                    PreferredDir[j]=M4DirCurr;
                    CurrF->SetS();

                    HasMoved=true;
                    break;
                }
            }
        }while (HasMoved);

        //then set the Faces finally
        IndexF.clear();
        for (size_t i=0;i<Mesh.face.size();i++)
        {
            if (!Mesh.face[i].IsS())continue;
            IndexF.push_back(i);
            //then make the direction on M2
            PreferredDir[i]=PreferredDir[i]%2;
            assert(PreferredDir[i]>=0);
            assert(PreferredDir[i]<4);
        }
    }

    //given the mesh and a set of sharp features pos return the array of faces on both sides
    static void GetFaceSharpPartitions(MeshType &Mesh,
                                       const std::vector<vcg::face::Pos<FaceType> > &SharpPos,
                                       std::vector<size_t> &IndexF0,
                                       std::vector<size_t> &IndexF1)
    {
        InitSharpFeaturesAsCreases(Mesh,SharpPos);

        vcg::face::Pos<FaceType> Pos0=SharpPos[0];

        //set the quaity of faces as -1
        vcg::tri::UpdateQuality<MeshType>::FaceConstant(Mesh,-1);

        //get sharp faces on first side
        GetFaceSharpPartitionSide(Mesh,SharpPos[0],IndexF0);

        //then go on the other side
        Pos0.FlipF();
        GetFaceSharpPartitionSide(Mesh,Pos0,IndexF1);
    }

    //given a pos return the closest principal direction
    static int getFaceEdgeOrientedDir(const vcg::face::Pos<FaceType> &SharpPos)
    {
        CoordType Dir=SharpPos.V()->P()-SharpPos.VFlip()->P();
        Dir.Normalize();

        //set preferred direction
        CoordType PD1=SharpPos.F()->PD1();
        CoordType PD2=SharpPos.F()->PD2();

        if (fabs(Dir*PD1)>fabs(Dir*PD2))
            return 0;
        else
            return 1;
    }

//    static void WhichM4Dir(const FaceType &F,
//                           const CoordType &Dir,
//                           int &PreferredDirM4)
//    {

//        CoordType PD[4];
//        //set preferred direction
//        PD[0]=f->PD1();
//        PD[1]=f->PD2();
//        PD[2]=-PD[0];
//        PD[3]=-PD[1];
//        ScalarType maxDot=-1;
//        PreferredDirM4=-1;
//        for (int i=0;i<4;i++)
//        {
//            if (Dir*PD[i]<maxDot)continue;
//            maxDot=Dir*PD[i];
//            PreferredDirM4=i;
//        }
//        assert(PreferredDirM4>=0);
//    }

//    static void getFaceEdgeM4DirNodes(AnisotropicGraph<MeshType> &AniGraph,
//                                      const FaceType &F,
//                                      const CoordType &Dir,
//                                      std::vector<size_t> &Nodes)
//    {
//        int PreferredDirM4;
//        int M4Dir=WhichM4Dir(F,Dir,PreferredDirM4);

//    }

    //given a mesh and a set of sharp features return the preferred dir
    static void getFaceCreasePreferredDir(MeshType &Mesh,
                                          const std::vector<vcg::face::Pos<FaceType> > &SharpPos,
                                          std::vector<int> &PreferredDir)
    {
        InitSharpFeaturesAsCreases(Mesh,SharpPos);
        PreferredDir=std::vector<int>(Mesh.face.size(),-1);
        for (size_t i=0;i< Mesh.face.size();i++)
            for (size_t j=0;j<3;j++)
            {
                //if (!Mesh.face[i].IsCrease(j))continue;
                if (!Mesh.face[i].IsFaceEdgeS(j))continue;
                FaceType *f=&Mesh.face[i];
                CoordType Dir=f->P0(j)-f->P1(j);
                Dir.Normalize();

                //set preferred direction
                CoordType PD1=f->PD1();
                CoordType PD2=f->PD2();

                if (fabs(Dir*PD1)>fabs(Dir*PD2))
                    PreferredDir[i]=0;
                else
                    PreferredDir[i]=1;
            }
    }

    //given a mesh and a set of sharp features return the faces on both sides and the preferred dir
    static void GetSingleSharpFacesInfo(MeshType &Mesh,
                                        const std::vector<vcg::face::Pos<FaceType> > &SharpPos,
                                        const std::vector<size_t> SharpEndP,
                                        std::vector<size_t> &IndexF0,
                                        std::vector<size_t> &IndexF1,
                                        std::vector<int> &PreferredDir)
    {
        //then set for each face that is on a sharp (concave) feature , the preferred direction
        //obviously there should be no face with two sharp edges
        getFaceCreasePreferredDir(Mesh,SharpPos,PreferredDir);

        //get the faces on both sides of the sharp feature
        GetFaceSharpPartitions(Mesh,SharpPos,IndexF0,IndexF1);

        //then propagate
        PropagatePropertiesOnSharpPartition(Mesh,SharpPos,SharpEndP,PreferredDir,IndexF0);
        PropagatePropertiesOnSharpPartition(Mesh,SharpPos,SharpEndP,PreferredDir,IndexF1);

        vcg::tri::UpdateFlags<MeshType>::VertexClearS(Mesh);
    }

    //get for each mesh and a set of sharp position return the set of
    static void GetSharpFacesInfo(MeshType &Mesh,
                                  const std::vector<std::vector<vcg::face::Pos<FaceType> > > &SharpPos,
                                  std::vector<std::vector<size_t> > &IndexF0,
                                  std::vector<std::vector<size_t> > &IndexF1,
                                  std::vector<int> &PreferredDir)
    {
        std::vector<size_t> EndPoints;


        SelectSharpEndPoints(Mesh,SharpPos);


        for (size_t i=0;i<Mesh.vert.size();i++)
            if (Mesh.vert[i].IsS())EndPoints.push_back(i);

        IndexF0.clear();
        IndexF1.clear();
        PreferredDir.clear();
        PreferredDir.resize(Mesh.face.size(),-1);
        for (size_t i=0;i<SharpPos.size();i++)
        {

            std::vector<size_t> IndexF0Temp,IndexF1Temp;
            std::vector<int> PreferredDirTemp;
            GetSingleSharpFacesInfo(Mesh,SharpPos[i],EndPoints,IndexF0Temp,IndexF1Temp,PreferredDirTemp);


            IndexF0.push_back(IndexF0Temp);
            IndexF1.push_back(IndexF1Temp);

            assert(PreferredDirTemp.size()==Mesh.face.size());
            for (size_t j=0;j<PreferredDirTemp.size();j++)
            {
                if (PreferredDirTemp[j]==-1)continue;
                if (PreferredDir[j]!=-1)
                {
                    vcg::tri::UpdateSelection<MeshType>::Clear(Mesh);
                    Mesh.face[j].SetS();
                    vcg::tri::io::ExporterPLY<MeshType>::Save(Mesh,"test.ply",vcg::tri::io::Mask::IOM_FACEFLAGS);
                    assert(0);
                }

                assert(PreferredDir[j]==-1);
                PreferredDir[j]=PreferredDirTemp[j];
            }
        }

    }

    //get end points of the sharp features
    static void GetSharpEndPoints(MeshType &Mesh,
                                  const std::vector<std::vector<vcg::face::Pos<FaceType> > > &SharpPos,
                                  std::vector<size_t> &SharpEndP)
    {
        SharpEndP.clear();
        for (size_t i=0;i<SharpPos.size();i++)
        {
            std::vector<size_t> SharpEndPTemp;
            GetSharpEndPoints(Mesh,SharpPos[i],SharpEndPTemp);
            SharpEndP.insert(SharpEndP.end(),SharpEndPTemp.begin(),SharpEndPTemp.end());
        }

        std::sort(SharpEndP.begin(),SharpEndP.end());
        std::vector<size_t>::iterator iteTerm=std::unique(SharpEndP.begin(),SharpEndP.end());
        SharpEndP.erase(iteTerm, SharpEndP.end());
    }

    //select edn points of the sharp features
    static void SelectSharpEndPoints(MeshType &Mesh,const std::vector<std::vector<vcg::face::Pos<FaceType> > > &SharpPos)
    {
        std::vector<size_t> SharpEndP;
        GetSharpEndPoints(Mesh,SharpPos,SharpEndP);
        for (size_t i=0;i<SharpEndP.size();i++)
            Mesh.vert[SharpEndP[i]].SetS();
    }

//    void SetCreases(MeshType &Mesh,const std::vector<vcg::face::Pos<FaceType> >  &SharpPos)
//    {
//        for (size_t j=0;j<SharpPos.size();j++)
//        {

//            vcg::face::Pos<FaceType> currPos=SharpPos[j];
//            size_t F0=vcg::tri::Index(Mesh,currPos.F());
//            size_t E0=currPos.E();
//            //go on the other side
//            currPos.FlipF();
//            size_t F1=vcg::tri::Index(Mesh,currPos.F());
//            size_t E1=currPos.E();

//            //set faux edges
//            Mesh.face[F0].SetCrease(E0);
//            Mesh.face[F1].SetCrease(E1);

//        }
//    }

//    void SetCreases(MeshType &Mesh,
//                    const std::vector<std::vector<vcg::face::Pos<FaceType> > > &SharpPos)
//    {
//        for (size_t i=0;i<SharpPos.size();i++)
//            SetCreases(Mesh,SharpPos[i]);
//    }

    //for each node returns the
    static void GetNodeConnectivityFactor(AnisotropicGraph<MeshType> &AniGraph,
                                          std::vector<size_t> &NodeConnectivityFactor)
    {
        std::vector<size_t> NodeIn(AniGraph.NumNodes(),0);
        std::vector<size_t> NodeOut(AniGraph.NumNodes(),0);
        for (size_t i = 0; i < AniGraph.Graph.size(); ++i)
            for (NeighborsIterator nIt = AniGraph.Graph[i]->neighbors.begin();
                 nIt != AniGraph.Graph[i]->neighbors.end(); ++nIt)
            {
                if (!nIt->Valid)continue;
                NodeOut[i]++;
                NodeIn[nIt->nodeID]++;
            }
        NodeConnectivityFactor.clear();
        NodeConnectivityFactor.resize(AniGraph.NumNodes(),0);
        for (size_t i=0;i<NodeOut.size();i++)
            NodeConnectivityFactor[i]=std::min(NodeOut[i],NodeIn[i]);
    }


    static int GetClosestLoop(AnisotropicGraph<MeshType> &anigraph,
                          std::vector<PathType> &Paths,
                          CoordType pos)
    {
        std::vector<bool> IsLoop(Paths.size(),true);
        std::vector<std::vector<std::pair<CoordType,CoordType> > > EdgePos;
        anigraph.GetPosPairs(Paths,IsLoop,EdgePos);
        int bestL=0;
        ScalarType minD=std::numeric_limits<ScalarType>::max();

        for (size_t i=0;i<EdgePos.size();i++)
            for (size_t j=0;j<EdgePos[i].size();j++)
            {
                vcg::Segment3<ScalarType> S(EdgePos[i][j].first,EdgePos[i][j].second);
                CoordType clos;
                ScalarType distTest;
                vcg::SegmentPointDistance(S,pos,clos,distTest);
                if (distTest>minD)continue;
                bestL=i;
                minD=distTest;
            }
        assert(bestL>=0);
        return bestL;
    }

    static void GetClosestLoopNode(AnisotropicGraph<MeshType> &anigraph,
                                  std::vector<PathType> &Paths,
                                  CoordType test_pos,int &IndexLoop,
                                  int &IndexNode)
    {
        ScalarType minD=std::numeric_limits<ScalarType>::max();
        IndexLoop=-1;IndexNode=-1;
        for (size_t i=0;i<Paths.size();i++)
            for (size_t j=0;j<Paths[i].nodes.size();j++)
            {
                CoordType pos=anigraph.NodePos(Paths[i].nodes[j]);
                ScalarType testDist=(pos-test_pos).Norm();
                if (testDist>minD)continue;
                minD=testDist;
                IndexLoop=i;
                IndexNode=j;
            }
        assert(IndexLoop>=0);
        assert(IndexNode>=0);
    }

    static ScalarType ComputeLoopCurvature(AnisotropicGraph<MeshType> &anigraph,
                                    const PathType &Path,
                                    const std::vector<CoordType> &CurvDir0,
                                    const std::vector<CoordType> &CurvDir1,
                                    const std::vector<ScalarType> &K0,
                                    const std::vector<ScalarType> &K1)
    {
         ScalarType CurvLoops=0;
         for (size_t i=0;i<Path.nodes.size();i++)
         {
             int IdxNode0=Path.nodes[i];
             int IdxNode1=Path.nodes[(i+1)%Path.nodes.size()];
             CoordType P0=anigraph.NodePos(IdxNode0);
             CoordType P1=anigraph.NodePos(IdxNode1);
             CoordType Dir=P1-P0;
             Dir.Normalize();
             size_t IndexF,IndexE,IndexV,M4Dir;
             anigraph.GetFaceEdgeDir(IdxNode0,IndexF,IndexE,M4Dir);
             IndexV=vcg::tri::Index(anigraph.Mesh(),anigraph.Mesh().face[IndexF].V(IndexE));

                if (fabs(Dir*CurvDir1[IndexV])>fabs(Dir*CurvDir0[IndexV]))
                    CurvLoops+=K0[IndexV];
                else
                    CurvLoops+=K1[IndexV];
            }
          CurvLoops/=Path.nodes.size();
          return CurvLoops;
    }

    static void CheckOppositeValid(AnisotropicGraph<MeshType> &anigraph)
    {
        for (size_t i=0;i<anigraph.NumNodes();i++)
        {
            if(!anigraph.IsValid(i))continue;
            int OppNode=anigraph.OppositeNode(i);
            assert(OppNode>=0);
            assert(OppNode<anigraph.NumNodes());
            assert(anigraph.IsValid(OppNode));
        }
    }
};



template < class MeshType >
class ConflictTable
{
    typedef typename MeshType::CoordType    CoordType;
    typedef typename MeshType::VertexType   VertexType;
    typedef typename MeshType::FaceType     FaceType;
    typedef typename MeshType::ScalarType   ScalarType;
    typedef typename AnisotropicGraph<MeshType>::NeighborsIterator NeighborsIterator;

    std::vector<std::vector<int> > CrossNodes;
    std::vector<std::vector<int> > TangentNodes;
    std::vector<std::vector<std::vector<int> > > FaceNodes;

    typedef std::pair<int,int> EdgeType;

    struct KeyHasher
    {
      std::size_t operator()(const EdgeType& k) const
      {
        using std::size_t;
        using std::hash;
        using std::string;

          return ((hash<int>()(k.first)
                ^ (hash<int>()(k.second) << 1)) >> 1);
      }
    };

    //std::unordered_map<EdgeType,std::vector<EdgeType>,KeyHasher > CrossEdges;
    //std::unordered_map<EdgeType,std::vector<EdgeType>,KeyHasher> TangentEdges;

    typedef typename AnisotropicGraph<MeshType>::Path PathType;


    bool Intersect(AnisotropicGraph<MeshType> &anigraph,
                   const EdgeType &E0,
                   const EdgeType &E1)
    {
        CoordType P0[2];
        CoordType P1[2];
        P0[0]=anigraph.NodePos(E0.first);
        P0[1]=anigraph.NodePos(E0.second);
        P1[0]=anigraph.NodePos(E1.first);
        P1[1]=anigraph.NodePos(E1.second);
        if (P0[0]==P1[0])return true;
        if (P0[0]==P1[1])return true;
        if (P0[1]==P1[0])return true;
        if (P0[1]==P1[1])return true;

        //check real intersection
        CoordType Dir0=(P0[0]-P1[0])^(P0[0]-P1[1]);
        CoordType Dir1=(P0[1]-P1[0])^(P0[1]-P1[1]);
//        Dir0.Normalize();
//        Dir1.Normalize();

        if ((Dir0*Dir1)<0)return true;

        return false;
    }

    void InitVertTable(AnisotropicGraph<MeshType> &anigraph)
    {
        CrossNodes.clear();
        TangentNodes.clear();
        CrossNodes.resize(anigraph.NumNodes());
        TangentNodes.resize(anigraph.NumNodes());
        for (size_t i=0;i<anigraph.NumNodes();i++)
        {
            LoopFunctions<MeshType>::GetTangentNodes(anigraph,i,TangentNodes[i]);
            LoopFunctions<MeshType>::GetCrossNodes(anigraph,i,CrossNodes[i]);
        }
    }

//    void InitEdgeTable(AnisotropicGraph<MeshType> &anigraph)
//    {
//        CrossEdges.clear();
//        TangentEdges.clear();

//        std::vector<std::vector<EdgeType> > FaceEdges;
//        FaceEdges.resize(anigraph.Mesh().face.size());

//        for (size_t i=0;i<anigraph.NumNodes();i++)
//        {
//            size_t IndexN0=i;
//            size_t IndexF,M4Dir;
//            anigraph.GetFaceDir(IndexN0,IndexF,M4Dir);

//            std::vector<size_t> NeighNodes,FacesId,M4Dirs;
//            std::vector<bool> Valid;
//            anigraph.GetNeighInfo(i,NeighNodes,FacesId,M4Dirs,Valid);
//            for (int j=0;j<NeighNodes.size();j++)
//            {
//                size_t IndexN1=NeighNodes[j];
//                EdgeType E(IndexN0,IndexN1);
//                FaceEdges[IndexF].push_back(E);
//            }
//        }

//        for (size_t i=0;i<FaceEdges.size();i++)
//        {
//            if (FaceEdges[i].size()<2)continue;
//            for (size_t j=0;j<FaceEdges[i].size()-1;j++)
//                for (size_t k=(j+1);k<FaceEdges[i].size();k++)
//                {
//                    EdgeType E0= FaceEdges[i][j];
//                    EdgeType E1= FaceEdges[i][k];
//                    size_t M4Dir0=anigraph.NodeDir(E0.first);
//                    size_t M4Dir1=anigraph.NodeDir(E1.first);
//                    bool IsTangent=((M4Dir0 % 2)==(M4Dir1 % 2));
//                    //if (!Intersect(anigraph,E0,E1))continue;
//                    if (IsTangent)
//                    {
//                        TangentEdges[E0].push_back(E1);
//                        TangentEdges[E1].push_back(E0);
//                    }
//                    else
//                    {
//                        CrossEdges[E0].push_back(E1);
//                        CrossEdges[E1].push_back(E0);
//                    }
//                }
//        }
//    }

    void InitFaceNodeTable(AnisotropicGraph<MeshType> &anigraph)
    {
        FaceNodes.clear();

        FaceNodes.resize(anigraph.Mesh().face.size());
        for (size_t i=0;i<FaceNodes.size();i++)
            FaceNodes[i].resize(2);

        for (size_t i=0;i<anigraph.NumNodes();i++)
        {
            size_t IndexF,M4Dir;
            anigraph.GetFaceDir(i,IndexF,M4Dir);

            FaceNodes[IndexF][M4Dir%2].push_back(i);
        }

    }

    void InvalidateTangantNodes(AnisotropicGraph<MeshType> &anigraph,PathType currP)
    {
        for (size_t i=0;i<currP.nodes.size();i++)
        {
            int IndexN=currP.nodes[i];
            std::vector<int> TangNodes;
            LoopFunctions<MeshType>::GetTangentNodes(anigraph,IndexN,TangNodes);
            for (size_t j=0;j<TangNodes.size();j++)
                anigraph.SetValid(TangNodes[j],false);
        }
    }

public:

    void InvalidateTangantEdges(AnisotropicGraph<MeshType> &anigraph,PathType currP,bool IsLoop=true)
    {
        int Limit=currP.nodes.size();
        if (!IsLoop)Limit--;
        for (int i=0;i<Limit;i++)
        {
            int IndexN0=currP.nodes[i];
            int IndexN1=currP.nodes[(i+1)%currP.nodes.size()];
            size_t IndexF,M4Dir;
            anigraph.GetFaceDir(IndexN0,IndexF,M4Dir);
            M4Dir=M4Dir%2;
            for (size_t  j=0;j<FaceNodes[IndexF][M4Dir].size();j++)
            {
                int IndexN2=FaceNodes[IndexF][M4Dir][j];
                if (!anigraph.IsValid(IndexN2))continue;

                //then get neighbors
                for (NeighborsIterator NeighIte=anigraph.Graph[IndexN2]->neighbors.begin();
                     NeighIte!=anigraph.Graph[IndexN2]->neighbors.end();NeighIte++)
                {
                    if ((*NeighIte).Valid==false)continue;
                    size_t IndexN3=(*NeighIte).nodeID;
                    EdgeType E0(IndexN0,IndexN1);
                    EdgeType E1(IndexN2,IndexN3);
                    if (!Intersect(anigraph,E0,E1))continue;
                    (*NeighIte).Valid=false;
                }
            }
        }
    }

    void GetConflictNodes(PathType currP,std::vector<int> &ConflNodes)
    {
        for (size_t i=0;i<currP.nodes.size();i++)
        {
            int IndexN=currP.nodes[i];
            ConflNodes.insert(ConflNodes.end(),TangentNodes[IndexN].begin(),TangentNodes[IndexN].end());
        }
    }


    void Init(AnisotropicGraph<MeshType> &anigraph)
    {
        int t0=clock();
        std::cout<<"Initializing Vert Intersection Table"<<std::endl;
        InitVertTable(anigraph);
        int t1=clock();
        std::cout<<"Time "<<t1-t0<<std::endl;
        std::cout<<"Initializing face node table"<<std::endl;
        InitFaceNodeTable(anigraph);
        std::cout<<"done"<<std::endl;
        int t2=clock();
        std::cout<<"Time "<<t2-t1<<std::endl;

    }

    void AddBarrierNonLoop(AnisotropicGraph<MeshType> &anigraph,
                            PathType &AddedPath)
    {

        InvalidateTangantNodes(anigraph,AddedPath);

        anigraph.InvalidateArcsOnNonValidNodes();

        InvalidateTangantEdges(anigraph,AddedPath,false);
    }

    void BarrierPaths(AnisotropicGraph<MeshType> &anigraph,
                      const std::vector<PathType> &TestLoops,
                      bool reset=true)
    {
        if (reset)
            anigraph.SetAllValid();

        for (size_t i=0;i<TestLoops.size();i++)
            InvalidateTangantNodes(anigraph,TestLoops[i]);

        anigraph.InvalidateArcsOnNonValidNodes();

        for (size_t i=0;i<TestLoops.size();i++)
            InvalidateTangantEdges(anigraph,TestLoops[i]);
    }


};



#endif //LOOP_OPTIMIZER
