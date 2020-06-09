#ifndef LOOP_OPTIMIZER
#define LOOP_OPTIMIZER

//#include <anisotropic_geodesic.h>
#include <tracing_field/graph/anisotropic_graph.h>
#include <tracing_field/conflict_finder.h>
//#include "scip_solver.h"
//#include <tracing_field/conflict_finder.h>
#include <vcg/complex/algorithms/point_sampling.h>
#include <vcg/complex/allocate.h>
#include <vcg/complex/algorithms/closest.h>
#include <vcg/space/index/grid_static_ptr.h>
#include <tracing_field/anisotropic_geodesic.h>
//#include <vcg/complex/algorithms/curve_on_manifold.h>
#include "tracing_field/remesh/edge_mesh_type.h"
#include "tracing_field/loop_common_functions.h"
#include "tracing_field/sharp_feature_sampler.h"
#include "tracing_field/sharp_feature_manager.h"
#include "tracing_field/inter_cross_validator.h"
//#include <tracing_field/loop_relax.h>

#define Intersection_needs 4
/**
 * @brief The class that performs parametrization by considering the connections within the graph
 */
template < class MeshType >
class LoopFinder {
public:

    LoopFinder(AnisotropicGraph<MeshType> &_anigraph,
               SharpFeaturesManager<MeshType> &_sharp) : anigraph(_anigraph),sharp(_sharp)
    {
        SharpSampl=NULL;
        selected=-1;
        selected_loop=-1;
        selected_node=-1;
    }

    ~LoopFinder()
    {
        delete(SharpSampl);
    }

private:
    typedef typename MeshType::CoordType    CoordType;
    typedef typename MeshType::VertexType   VertexType;
    typedef typename MeshType::FaceType     FaceType;
    typedef typename MeshType::ScalarType   ScalarType;

    AnisotropicGraph<MeshType> &anigraph;               //the mesh on which the separatrix graph is calculated
    SharpFeaturesManager<MeshType> &sharp;

    typedef typename AnisotropicGraph<MeshType>::Path PathType;
    typedef typename AnisotropicGraph<MeshType>::NodeListIterator NodeListIterator;
    typedef typename AnisotropicGraph<MeshType>::NodeListConstIterator NodeListConstIterator;

    typedef typename ConflictFinder<MeshType>::ConflictInfo ConflictInfo;

    //they are removed because of kernel distance
    std::vector<bool> ActiveNodes;
    std::vector<size_t> SingNodes;
    std::set<size_t> SingNodesSet;

    bool AddingEssential;

    enum LoopType{GuardLoop,RegularLoop,ConcaveLoop};

    int selected;

public:
    int selected_loop;
    int selected_node;

private:
    struct CandidateLoop
    {
        //if need update or not
        bool To_Update;

        //if the initial nod or other has been already traced
        bool BunrdOut;

        //initial tracing node
        int InitialNode;

        //the kind of loop
        PathType LoopPath;

        //the averagw distance wrt other loops
        ScalarType AvgDistance;

        //the list of sharp features and loops it intersects
        std::vector<int> SharpInters,LoopInters;

        //the intersections positions
        std::vector<CoordType> SharpIntersPos,LoopIntersPos;

        //number of sharp features that the candidate solve
        int SolveSharp;

        //the number of loops that the candidate solve
        int SolveLoop;

        //the number of inter cross loops that the candidate solve
        int SolveInterCross;

        //the self distance
        ScalarType SelfDist;

        //the generative nodes used when tracing
        std::vector<size_t> OrthoGenerativeNodes;

        //the tangent nodes used as a guard or to compute distance
        std::vector<std::vector<int> > TangentGuardNodes;

        //TODO REMOVE THAT
        LoopType LType;
        bool IsLoop;

        //if is essential or not
        bool IsEssential;

        CandidateLoop()
        {
            InitialNode=-1;
            AvgDistance=-1;
            LType=RegularLoop;
            IsLoop=true;
            IsEssential=false;
            SolveSharp=0;
            SolveLoop=0;
            SolveInterCross=0;
        }

        CandidateLoop(int _InitialNode)
        {
            InitialNode=_InitialNode;
            AvgDistance=-1;
            BunrdOut=false;
            To_Update=true;
            LType=RegularLoop;
            IsLoop=true;
            SolveSharp=0;
            SolveLoop=0;
            SolveInterCross=0;
            IsEssential=false;
            //OpenNum=0;
            //SuperfluousNum=0;

            //SelfDist=0;
        }

        bool IsValid()const
        {
            if (AvgDistance<0){std::cout<<"worng avg distance "<<AvgDistance<<std::endl;return false;}
            if (InitialNode<0){std::cout<<"worng init node "<<InitialNode<<std::endl;return false;}
            if (BunrdOut){std::cout<<"burn out init node "<<std::endl;return false;}
            return true;
        }

        inline bool operator <(const CandidateLoop &Cloop)const
        {
            assert(Cloop.IsValid());
            assert(IsValid());

            bool SolveSh0=(SolveSharp>0);
            bool SolveSh1=(Cloop.SolveSharp>0);
            bool SolveL0=(SolveLoop>0);
            bool SolveL1=(Cloop.SolveLoop>0);

            if (SolveSh0 && (!SolveSh1))return false;
            if (SolveSh1 && (!SolveSh0))return true;
            if (SolveSh0 && SolveSh1)return AvgDistance<Cloop.AvgDistance;

            if (SolveL0 && (!SolveL1))return false;
            if (SolveL1 && (!SolveL0))return true;

            if (SolveInterCross!=Cloop.SolveInterCross)
                return (SolveInterCross<Cloop.SolveInterCross);

            return AvgDistance<Cloop.AvgDistance;
        }

        //                inline bool operator <(const CandidateLoop &Cloop)const
        //                {
        //                    assert(Cloop.IsValid());
        //                    assert(IsValid());

        //                    if (SolveSharp !=Cloop.SolveSharp)return (SolveSharp <Cloop.SolveSharp);
        //                    if (SolveLoop !=Cloop.SolveLoop)return (SolveLoop <Cloop.SolveLoop);

        //                    return AvgDistance<Cloop.AvgDistance;
        //                }

        inline void ClearNodes()
        {
            LoopPath.clear();
        }

        inline void RestoreNodes(std::vector<size_t> &NodeIndexes)
        {
            LoopPath.nodes=NodeIndexes;
        }
    };

    std::vector<CandidateLoop> Candidates;
    std::vector<CandidateLoop> ChoosenLoops;
    //std::vector<bool> SelectedLoops;
    std::vector<PathType> ConvexPath;

    struct SharpFeature
    {
        //the set of pos defining the sharp fearure
        std::vector<vcg::face::Pos<FaceType> > FeaturePos;

        //if is loop or nor
        bool IsLoop;

        //concave or not
        bool IsConcave;

        //the set of guard nodes
        std::vector<size_t> TangentGuardNodes;

        //the tangent nodes used as a guard
        std::vector<size_t> OrthoGenerativeNodes;

        //the set of positions used to check if is intersected
        std::set<CoordType> IntersectingPos;

        SharpFeature(const std::vector<vcg::face::Pos<FaceType> > &_FeaturePos,
                     bool _IsLoop,bool _IsConcave)
        {
            IsLoop=_IsLoop;
            IsConcave=_IsConcave;
            FeaturePos=std::vector<vcg::face::Pos<FaceType> >(_FeaturePos.begin(),_FeaturePos.end());
        }
    };

    std::vector<SharpFeature> SharpFeatures;
    std::vector<ScalarType> NodeDist;

    size_t InitialLoopNum;
    size_t UpdatedLoopNum;

    size_t TimeFirstLoops;
    size_t TimeFirstGuess;
    size_t TimeLoopUpdate;
    size_t TimePriorityUpdate;

public:

    void GetFinalPaths(std::vector<PathType> &ChoosenPaths,
                       std::vector<bool> &IsLoop,
                       std::vector<bool> &Essential,
                       std::vector<bool> &CrossOK)
    {
        //Essential
        ChoosenPaths.clear();
        IsLoop.clear();
        Essential.clear();
        CrossOK.clear();
        for (size_t i=0;i<ChoosenLoops.size();i++)
        {
            if (ChoosenLoops[i].LType==GuardLoop)continue;
            ChoosenPaths.push_back(ChoosenLoops[i].LoopPath);
            IsLoop.push_back(ChoosenLoops[i].IsLoop);
            Essential.push_back(ChoosenLoops[i].IsEssential);

            if (LoopIntersectionNeeds[i]==0)
                CrossOK.push_back(true);
            else
                CrossOK.push_back(false);
        }

        for (size_t i=0;i<ConvexPath.size();i++)
        {
            ChoosenPaths.push_back(ConvexPath[i]);
            IsLoop.push_back(false);
            Essential.push_back(false);
            CrossOK.push_back(true);
        }

    }

    void GetConvexPath(std::vector<PathType> &_ConvexPath)
    {
        _ConvexPath=ConvexPath;
    }

    void GetChoosenPaths(std::vector<PathType> &ChoosenPaths)
    {
        ChoosenPaths.clear();
        for (size_t i=0;i<ChoosenLoops.size();i++)
            ChoosenPaths.push_back(ChoosenLoops[i].LoopPath);
    }



private:

    MeshType VertMesh;
    vcg::GridStaticPtr<VertexType,ScalarType> GridVert;

    void UpdatePriorityBySelfDist(CandidateLoop &CLoop)
    {
        assert(!CLoop.BunrdOut);
        assert(CLoop.LoopPath.nodes.size()>0);

        std::vector<PathType> TestLoop;
        TestLoop.push_back(CLoop.LoopPath);
        ScalarType SelfD=AnisotropicQuery<MeshType>::LoopSelfTangentialDist(anigraph,TestLoop,SParams.MaxDev,SParams.PenaltyDrift,2);
        CLoop.AvgDistance=SelfD;
    }

    void UpdatePriorityBySelfDist()
    {
        for (size_t i=0;i<Candidates.size();i++)
            UpdatePriorityBySelfDist(Candidates[i]);
    }
    //    void UpdatePriorityByLoopDistances(CandidateLoop &CLoop,
    //                                       const std::vector<ScalarType> &NodeDist,
    //                                       ScalarType MaxD=std::numeric_limits<ScalarType>::max())
    //    {
    //        assert(!CLoop.BunrdOut);
    //        assert(CLoop.LoopPath.nodes.size()>0);

    //        assert(NodeDist.size()==anigraph.NumNodes());
    //        ScalarType AvgDist=0;
    //        //ScalarType AvgDist=MaxD;
    //        for (size_t j=0;j<CLoop.LoopPath.nodes.size();j++)
    //        {
    //            size_t NodeI=CLoop.LoopPath.nodes[j];

    //            //            if (SingNodesSet.count(NodeI)>0)
    //            //            {
    //            //                CLoop.Priority=0;
    //            //                return;
    //            //            }
    //            assert(NodeI<NodeDist.size());
    //            AvgDist+=std::min(MaxD,NodeDist[NodeI]);
    //            //AvgDist=std::min(MaxD,NodeDist[NodeI]);
    //        }

    //        AvgDist/=(ScalarType)CLoop.LoopPath.nodes.size();

    //        CLoop.AvgDistance=AvgDist;//+CLoop.SelfDist);
    //    }

    void UpdatePriorityByLoopDistances(CandidateLoop &CLoop,
                                       const std::vector<ScalarType> &NodeDist)
    {
        assert(!CLoop.BunrdOut);
        assert(CLoop.LoopPath.nodes.size()>0);

        assert(NodeDist.size()==anigraph.NumNodes());
        ScalarType AvgDist=0;
        //ScalarType MinDist=std::numeric_limits<ScalarType>::max();

        for (size_t i=0;i<CLoop.TangentGuardNodes.size();i++)
        {
            assert(CLoop.TangentGuardNodes[i].size()>0);
            ScalarType MindD=std::numeric_limits<ScalarType>::max();
            for (size_t j=0;j<CLoop.TangentGuardNodes[i].size();j++)
            {
                int NodeI=CLoop.TangentGuardNodes[i][j];
                assert(NodeI>=0);
                assert(NodeI<(int)NodeDist.size());
                MindD=std::min(MindD,NodeDist[NodeI]);
            }
            AvgDist+=MindD;
            //MinDist=std::min(MinDist,MinD);
        }

        AvgDist/=(ScalarType)CLoop.LoopPath.nodes.size();

        CLoop.AvgDistance=AvgDist;//+CLoop.SelfDist);
        //CLoop.AvgDistance=MinDist;
    }

    void UpdateCandidatePriorityByLoopDistances()
    {
        size_t t3=clock();
        for (size_t i=0;i<Candidates.size();i++)
            UpdatePriorityByLoopDistances(Candidates[i],NodeDist);
        size_t t4=clock();
        TimePriorityUpdate+=t4-t3;
    }


    //    void UpdateCandidatePriorityByIntersections()
    //    {
    //        size_t t3=clock();
    //        for (size_t i=0;i<Candidates.size();i++)
    //        {
    //            std::vector<int> TestSharpIntersectionNeeds=SharpIntersectionNeeds;
    //            std::vector<int> TestLoopIntersectionNeeds=LoopIntersectionNeeds;

    //            std::vector<int> SharpInters,LoopInters;
    //            std::vector<CoordType> SharpIntersPos,LoopIntersPos;
    //            GetCandidateIntersections(Candidates[i],SharpInters,SharpIntersPos,LoopInters,LoopIntersPos);

    //            Candidates[i].SolveSharp=0;
    //            Candidates[i].SolveLoop=0;
    //            Candidates[i].OpenNum=4;
    //            Candidates[i].SuperfluousNum=0;

    //            //compute how much it solve
    //            for (size_t j=0;j<SharpInters.size();j++)
    //            {
    //                int IndexShF=SharpInters[j];
    //                assert(IndexShF<TestSharpIntersectionNeeds.size());
    //                assert(IndexShF>=0);

    //                if (SharpFeatures[IndexShF].IsConcave)Candidates[i].OpenNum--;

    //                if (TestSharpIntersectionNeeds[IndexShF]>0)
    //                {
    //                    TestSharpIntersectionNeeds[IndexShF]--;
    //                    Candidates[i].SolveSharp++;
    //                }
    //                else
    //                {
    //                    if (SharpFeatures[IndexShF].IsLoop)
    //                        Candidates[i].SuperfluousNum++;
    //                }
    //            }

    //            for (size_t j=0;j<LoopInters.size();j++)
    //            {
    //                Candidates[i].OpenNum--;
    //                int IndexL=LoopInters[j];
    //                assert(IndexL<TestLoopIntersectionNeeds.size());
    //                assert(IndexL>=0);
    //                if (TestLoopIntersectionNeeds[IndexL]>0)
    //                {
    //                    TestLoopIntersectionNeeds[IndexL]--;
    //                    Candidates[i].SolveLoop++;
    //                }
    //                //                else
    //                //                    Candidates[i].SuperfluousNum++;
    //            }

    //            //and then how many it opens
    //            //Candidates[i].OpenNum= 4- LoopInters.size() - SharpInters.size();
    //            Candidates[i].OpenNum= 4- std::max(Candidates[i].OpenNum,0);

    //            //then  save the one that are affected
    //            Candidates[i].SharpInters=SharpInters;
    //            Candidates[i].LoopInters=LoopInters;
    //            Candidates[i].SharpIntersPos=SharpIntersPos;
    //            Candidates[i].LoopIntersPos=LoopIntersPos;
    //        }

    //        size_t t4=clock();
    //        TimePriorityUpdate+=t4-t3;
    //    }

    void UpdateCandidatePriorityByIntersections()
    {
        size_t t3=clock();
        for (size_t i=0;i<Candidates.size();i++)
        {
            //update if needed
            if (Candidates[i].To_Update)
            {
                std::vector<int> SharpInters,LoopInters;
                std::vector<CoordType> SharpIntersPos,LoopIntersPos;
                GetCandidateIntersections(Candidates[i],SharpInters,SharpIntersPos,LoopInters,LoopIntersPos);
                //then  save the one that are affected
                Candidates[i].SharpInters=SharpInters;
                Candidates[i].LoopInters=LoopInters;
                Candidates[i].SharpIntersPos=SharpIntersPos;
                Candidates[i].LoopIntersPos=LoopIntersPos;
            }

            //then compute
            std::vector<int> TestSharpIntersectionNeeds=SharpIntersectionNeeds;
            std::vector<int> TestLoopIntersectionNeeds=LoopIntersectionNeeds;


            Candidates[i].SolveSharp=0;
            Candidates[i].SolveLoop=0;
            //Candidates[i].OpenNum=4;
            //Candidates[i].SuperfluousNum=0;

            //compute how much it solve
            for (size_t j=0;j<Candidates[i].SharpInters.size();j++)
            {
                int IndexShF=Candidates[i].SharpInters[j];
                assert(IndexShF<(int)TestSharpIntersectionNeeds.size());
                assert(IndexShF>=0);

                //if (SharpFeatures[IndexShF].IsConcave)Candidates[i].OpenNum--;

                if (TestSharpIntersectionNeeds[IndexShF]>0)
                {
                    TestSharpIntersectionNeeds[IndexShF]--;
                    Candidates[i].SolveSharp++;
                }/*
                else
                {
                    //if (SharpFeatures[IndexShF].IsLoop)
                    Candidates[i].SuperfluousNum++;
                }*/
            }

            for (size_t j=0;j<Candidates[i].LoopInters.size();j++)
            {
                //Candidates[i].OpenNum--;
                int IndexL=Candidates[i].LoopInters[j];
                assert(IndexL<(int)TestLoopIntersectionNeeds.size());
                assert(IndexL>=0);
                if (TestLoopIntersectionNeeds[IndexL]>0)
                {
                    TestLoopIntersectionNeeds[IndexL]--;
                    Candidates[i].SolveLoop++;
                }
            }
        }

        size_t t4=clock();
        TimePriorityUpdate+=t4-t3;
    }

    void ComputeLoop(CandidateLoop &CLoop)
    {
        assert(!CLoop.BunrdOut);

        if (!ActiveNodes[CLoop.InitialNode])
        {
            CLoop.BunrdOut=true;
            return;
        }

        assert(anigraph.Graph[CLoop.InitialNode]->Valid);
        bool extracted=AnisotropicQuery<MeshType>::ExtractLoop(anigraph,CLoop.InitialNode,SParams.MaxDev,SParams.PenaltyDrift,CLoop.LoopPath);
        if (!extracted)
        {
            CLoop.BunrdOut=true;
            return;
        }
        ConflictFinder<MeshType> ConflFinder(anigraph);
        bool self_conflict=ConflFinder.SelfConflict(CLoop.LoopPath,2);
        if (self_conflict)
        {
            CLoop.BunrdOut=true;
            return;
        }

        if (SParams.AvoidSelfIntersections)
        {
            bool cross_int=ConflFinder.CrossIntersect(CLoop.LoopPath,true);
            if (cross_int)
            {
                CLoop.BunrdOut=true;
                return;
            }
        }
    }

    //    void UpdateCandidatesSelfDistances()
    //    {
    //        anigraph.SetAllValid();
    //        for (size_t i=0;i<Candidates.size();i++)
    //        {
    //            if (!Candidates[i].To_Update)continue;
    //            std::vector<PathType> TestLoop;
    //            TestLoop.push_back(Candidates[i].LoopPath);
    //            Candidates[i].SelfDist=AnisotropicQuery<MeshType>::(anigraph,TestLoop,45,0,2);
    //        }
    //    }

    void SetUpdateFlag(CandidateLoop &CLoop)
    {
        assert(CLoop.InitialNode>=0);
        CLoop.To_Update=false;

        //in this case is the first exectution
        if (CLoop.AvgDistance<0)
        {
            assert(CLoop.LoopPath.nodes.size()==0);
            CLoop.To_Update=true;
            return;
        }

        //otherwise check if need to be updated
        //        for (size_t i=0;i<CLoop.LoopPath.nodes.size();i++)
        //        {

        //            int NodeI=CLoop.LoopPath.nodes[i];
        //            assert(NodeI>=0);
        //            assert(NodeI<(int)ActiveNodes.size());
        //            if (ActiveNodes[NodeI])continue;
        //            if (SParams.UpdateLoops)
        //                CLoop.To_Update=true;
        //            else
        //                CLoop.BunrdOut=true;
        //            return;
        //        }
        for (size_t i=0;i<CLoop.TangentGuardNodes.size();i++)
        {
            for (size_t j=0;j<CLoop.TangentGuardNodes[i].size();j++)
            {
                int NodeI=CLoop.TangentGuardNodes[i][j];
                assert(NodeI>=0);
                assert(NodeI<(int)ActiveNodes.size());
                if (ActiveNodes[NodeI])continue;
                if (SParams.UpdateLoops)
                    CLoop.To_Update=true;
                else
                    CLoop.BunrdOut=true;
                return;
            }
        }
    }

    void SetBarrierNonActiveNodes()
    {

        anigraph.SetAllValid();
        assert(ActiveNodes.size()==anigraph.NumNodes());
        for (size_t i=0;i<ActiveNodes.size();i++)
        {
            if (ActiveNodes[i])continue;
            //get all nodes with the coordinate
            anigraph.SetValid(i,false);
        }

        //then disable all arcs that goes or exit from non valid nodes
        anigraph.InvalidateArcsOnNonValidNodes();

    }

    //    void SetBarrierChoosenLoops(bool Reset=true)
    //    {
    //        if (Reset)
    //            ActiveNodes.assign(ActiveNodes.size(),true);

    //        std::vector<PathType> ChosenP;
    //        GetChoosenPaths(ChosenP);
    //        assert(ChosenP.size()==ChoosenLoops.size());

    //        std::set<std::pair<size_t,size_t> > FaceDir;
    //        for (size_t i=0;i<ChosenP.size();i++)
    //        {
    //            if (ChoosenLoops[i].LType==GuardLoop)
    //            {
    //                for (size_t j=0;j<ChosenP[i].nodes.size();j++)
    //                {
    //                    int IndexNode=ChosenP[i].nodes[j];
    //                    assert(anigraph.IsEdgeNode(IndexNode));
    //                    size_t FaceId;
    //                    size_t M4Dir;
    //                    anigraph.GetFaceDir(IndexNode,FaceId,M4Dir);
    //                    FaceDir.insert(std::pair<size_t,size_t>(FaceId,M4Dir));
    //                    FaceDir.insert(std::pair<size_t,size_t>(FaceId,(M4Dir+2)%4));
    //                }
    //            }
    //            else
    //            {
    //                std::vector<std::pair<size_t,size_t> > CurrFaceDir;
    //                anigraph.GetFacesDir(ChosenP[i],CurrFaceDir,ChoosenLoops[i].IsLoop);
    //                for (size_t j=0;j<CurrFaceDir.size();j++)
    //                {
    //                    int FaceId=CurrFaceDir[j].first;
    //                    int M4Dir=CurrFaceDir[j].second;
    //                    FaceDir.insert(std::pair<size_t,size_t>(FaceId,M4Dir));
    //                    FaceDir.insert(std::pair<size_t,size_t>(FaceId,(M4Dir+2)%4));
    //                }
    //            }
    //        }

    //        //then invalidate all nodes considering such directions
    //        for (size_t i=0;i<anigraph.NumNodes();i++)
    //        {
    //            size_t FaceId;
    //            size_t M4Dir;
    //            if (!ActiveNodes[i])continue;
    //            if (!anigraph.IsEdgeNode(i))continue;
    //            anigraph.GetFaceDir(i,FaceId,M4Dir);
    //            std::pair<size_t,size_t> key(FaceId,M4Dir);
    //            if (FaceDir.count(key)==0)continue;
    //            ActiveNodes[i]=false;

    //            int IndexOpp=anigraph.OppositeNode(i);
    //            if ((SParams.HasBorder)&&(IndexOpp==-1))continue;
    //            assert(IndexOpp>=0);
    //            assert(IndexOpp<anigraph.NumNodes());
    //            assert(anigraph.IsEdgeNode(IndexOpp));
    //            ActiveNodes[IndexOpp]=false;
    //        }
    //        //        if (SParams.StayAwayFromSing)
    //        //        {
    //        //            for (size_t i=0;i<SingGuardNodes.size();i++)
    //        //            {
    //        //                int IndexN=SingGuardNodes[i];
    //        //                assert(anigraph.IsEdgeNode(IndexN));
    //        //                ActiveNodes[IndexN]=false;

    //        //            }
    //        //        }
    //        if (Reset)
    //            SetBarrierNonActiveNodes();
    //    }

    void EraseBurndOutCandidates()
    {
        //erase the bourndout
        std::vector<CandidateLoop> TempCandidates;
        for (size_t i=0;i<Candidates.size();i++)
        {
            if (Candidates[i].BunrdOut)continue;
            TempCandidates.push_back(Candidates[i]);
        }
        Candidates=std::vector<CandidateLoop>(TempCandidates.begin(),TempCandidates.end());
    }

    //    void UpdateStep()
    //    {
    //        //first get current configuration and uptate distances and valid
    //        //std::vector<ScalarType> NodeDist;
    //        anigraph.SetAllValid();

    //        std::cout<<"Updating step"<<std::endl;

    //        //if there are already loops compute distances from them
    //        if (ChoosenLoops.size()>0)
    //        {

    //            std::vector<PathType> ChoosenPaths;
    //            GetChoosenPaths(ChoosenPaths);
    //            size_t t0=clock();

    //            std::cout<<"Computing loops distances"<<std::endl;
    //            AnisotropicQuery<MeshType>::ComputeLoopDistances(anigraph,ChoosenPaths,45,0,NodeDist);

    //            size_t t1=clock();
    //            TimePriorityUpdate+=t1-t0;

    //            std::cout<<"Setting non active nodes"<<std::endl;

    //            SetBarrierChoosenLoops();
    //        }

    //        //then set the ones to be updated
    //        int to_be_updated=0;
    //        for (size_t i=0;i<Candidates.size();i++)
    //        {
    //            SetUpdateFlag(Candidates[i]);
    //            if (Candidates[i].To_Update)
    //                to_be_updated++;
    //        }

    //        std::cout<<" Loops need to be updated "<< to_be_updated<<" on total of "<<Candidates.size()<<std::endl;

    //        //then update the loops
    //        size_t t0=clock();
    //        std::cout<<"*** UPDATING LOOPS ***"<<std::endl;
    //        size_t UpdatedNum=0;
    //        for (size_t i=0;i<Candidates.size();i++)
    //        {
    //            if (!Candidates[i].To_Update)continue;
    //            ComputeLoop(Candidates[i]);
    //            UpdatedNum++;
    //        }


    //        std::cout<<"*** End Updating LOOPS ***"<<std::endl;

    //        size_t t1=clock();
    //        if (ChoosenLoops.size()==0)
    //        {
    //            TimeFirstLoops+=t1-t0;
    //            InitialLoopNum=UpdatedNum;
    //        }
    //        else
    //        {
    //            TimeLoopUpdate+=t1-t0;
    //            UpdatedLoopNum+=UpdatedNum;
    //        }
    //        //trim the boundout
    //        EraseBurndOutCandidates();


    //        std::cout<<"There are "<<Candidates.size()<<" candidates"<<std::endl;
    //        //update priorities
    //        //std::cout<<"*** UPDATING LOOPS PRIORITIES ***"<<std::endl;
    //        anigraph.SetAllValid();

    //        if (ChoosenLoops.size()>0)
    //        {
    //            size_t t3=clock();
    //            for (size_t i=0;i<Candidates.size();i++)
    //                UpdatePriorityByLoopDistances(Candidates[i],NodeDist);
    //            size_t t4=clock();
    //            TimePriorityUpdate+=t4-t3;
    //        }
    //        else
    //        {
    //            int MaxEsteem=50;
    //            int Interval=floor(0.5+Candidates.size()/(ScalarType)MaxEsteem);
    //            Interval=std::max(1,Interval);
    //            size_t t3=clock();
    //            for (size_t i=0;i<Candidates.size();i++)
    //            {
    //                if ((i%Interval==0))
    //                    UpdatePriorityBySelfDist(Candidates[i]);
    //                else
    //                    Candidates[i].Priority=0.00000001;
    //            }
    //            size_t t4=clock();
    //            TimeFirstGuess+=t4-t3;
    //        }

    //    }

    void UpdateCandidates()
    {
        std::cout<<"Updating step"<<std::endl;

        //then update the loops
        //        size_t t0=clock();
        //anigraph.WeightLenghtByVertexQ();

        std::cout<<"*** UPDATING LOOPS ***"<<std::endl;
        size_t UpdatedNum=0;
        for (size_t i=0;i<Candidates.size();i++)
        {
            if (!Candidates[i].To_Update)continue;
            ComputeLoop(Candidates[i]);
            //increase number of updated
            UpdatedNum++;
        }
        std::cout<<" updated "<< UpdatedNum<<" on total of "<<Candidates.size()<<std::endl;
        std::cout<<"*** End Updating LOOPS ***"<<std::endl;

        //trim the boundout
        EraseBurndOutCandidates();

        //compute the guard nodes
        for (size_t i=0;i<Candidates.size();i++)
        {
            if (!Candidates[i].To_Update)continue;
            LoopFunctions<MeshType>::GetTangentNodes(anigraph,Candidates[i].LoopPath,Candidates[i].TangentGuardNodes);

            //increase number of updated
            UpdatedNum++;
        }


        std::cout<<"There are "<<Candidates.size()<<" candidates"<<std::endl;
    }

    bool SampleCandidates(size_t NumSeeds)
    {

        std::cout << "Sampling  Initial Points "<< std::endl<< std::flush;

        //first sample a set of evenly distributed nodes
        unsigned int randSeed=1232;
        std::vector<CoordType> poissonSamples;
        ScalarType Radius=vcg::tri::ComputePoissonDiskRadius(anigraph.Mesh(),NumSeeds);
        vcg::tri::PoissonSampling<MeshType>(anigraph.Mesh(),poissonSamples,NumSeeds,Radius,1,0.04f,randSeed);

        //then transform graph nodes to a mesh
        //        MeshType VertMesh;

        //vcg::tri::Allocator<MeshType>::AddVertices(VertMesh,anigraph.NumNodes());
        for (size_t i=0;i<anigraph.NumNodes();i++)
        {
            if (!ActiveNodes[i])continue;
            if (anigraph.IsBorderNode(i))continue;
            vcg::tri::Allocator<MeshType>::AddVertices(VertMesh,1);
            VertMesh.vert.back().P()=anigraph.NodePos(i);
        }

        vcg::tri::UpdateBounding<MeshType>::Box(VertMesh);
        //create a grid
        GridVert.Set(VertMesh.vert.begin(),VertMesh.vert.end());

        std::vector<size_t> seedNodes;
        //then for each sample retrieve the closeset node
        for (size_t i=0;i<poissonSamples.size();i++)
        {
            CoordType pos=poissonSamples[i];
            ScalarType minDist;
            VertexType *closV=vcg::tri::GetClosestVertex(VertMesh,GridVert,pos,std::numeric_limits<ScalarType>::max(),minDist);
            assert(closV!=NULL);

            //find the vertices and the directions relative
            CoordType Vpos=closV->P();

            std::vector<size_t> dirNodes;
            anigraph.GetPossibleDirNodes(Vpos,dirNodes,2);

            //add only nodes that are not deactivated
            for (size_t j=0;j<dirNodes.size();j++)
            {
                size_t currN=dirNodes[j];

                //this is for the next sampling step
                if (!ActiveNodes[currN])continue;

                //not starting a loop from singularities
                if (!anigraph.IsEdgeNode(currN))continue;

                assert(anigraph.IsValid(currN));

                seedNodes.push_back(currN);
            }
        }
        //sort and make initial nodes unique
        std::sort(seedNodes.begin(),seedNodes.end());

        seedNodes.erase( std::unique( seedNodes.begin(), seedNodes.end() ), seedNodes.end() );

        std::cout << "Sampled " << seedNodes.size() << " Initial Points "<< std::endl<< std::flush;
        //        //std::cout << "Already "<< DeactivatedNodes.size() << " inserted Nodes"<<std::endl;

        Candidates.clear();
        for (size_t i=0;i<seedNodes.size();i++)
            Candidates.push_back(CandidateLoop(seedNodes[i]));

        return (seedNodes.size()>0);
    }




    //std::vector<std::vector<vcg::face::Pos<FaceType> > > FeaturePos;
    //std::vector<bool> Concave;
    //std::vector<std::vector<size_t> > GuardNodes;

    //    void InitSharpFeaturesSimple()
    //    {

    //        FeaturePos.clear();

    //        //first find concave edges
    //        std::vector<std::vector<vcg::face::Pos<FaceType> > > ConcaveEdges;
    //        //anigraph.Mesh().FindSharpFeatures(ConcaveEdges,ConcaveEdge,SParams.SharpDegree);
    //        anigraph.Mesh().GetSharpFeatures(ConcaveEdges,ConcaveEdge);
    //        std::cout<<"* There are "<< ConcaveEdges.size()<<" convex sharp features" << std::endl;

    //        //then find convex edges
    //        std::vector<std::vector<vcg::face::Pos<FaceType> > > ConvexEdges;
    //        //anigraph.Mesh().FindSharpFeatures(ConvexEdges,ConvexEdge,SParams.SharpDegree);
    //        anigraph.Mesh().GetSharpFeatures(ConcaveEdges,ConvexEdge);
    //        std::cout<<"* There are "<< ConvexEdges.size()<<" convex sharp features" << std::endl;

    //        // then copy onto feature pos
    //        FeaturePos=ConcaveEdges;
    //        Concave.resize(FeaturePos.size(),true);

    //        FeaturePos.insert(FeaturePos.end(),ConvexEdges.begin(),ConvexEdges.end());
    //        Concave.resize(FeaturePos.size(),false);

    //        //NumSharpFeature=FeaturePos.size();

    //        std::cout<<"* Initializing portals " << std::endl;

    //    }

    void PrintFeaturesStats()
    {
        size_t NumConcave=0;
        size_t LoopConcave=0;
        size_t LoopConvex=0;
        for (size_t i=0;i<SharpFeatures.size();i++)
        {
            if (SharpFeatures[i].IsConcave)
            {
                NumConcave++;
                if (SharpFeatures[i].IsLoop)LoopConcave++;
            }else
            {
                if (SharpFeatures[i].IsLoop)LoopConvex++;
            }
        }
        std::cout<<"Total Features "<<SharpFeatures.size()<<std::endl;

        std::cout<<"Num Concave "<<NumConcave<<std::endl;
        std::cout<<"Num Concave Loops "<<LoopConcave<<std::endl;

        std::cout<<"Num Convex "<<SharpFeatures.size()-NumConcave<<std::endl;
        std::cout<<"Num Convex Loops "<<LoopConvex<<std::endl;
    }

    void SubsampleIntervals(std::vector<size_t> &Nodes,
                            size_t SampleDensity=20)
    {
        std::vector<size_t> UnifiedOrthoNode;

        //erase duplicates in position
        std::set<CoordType> NodePos;
        for (size_t i=0;i<Nodes.size();i++)
        {
            CoordType currP=anigraph.NodePos(Nodes[i]);
            if (NodePos.count(currP)>0)continue;
            NodePos.insert(currP);
            UnifiedOrthoNode.push_back(Nodes[i]);
        }

        if (UnifiedOrthoNode.size()<SampleDensity)
        {
            Nodes=UnifiedOrthoNode;
            return;
        }
        Nodes.clear();

        //then subsample
        ScalarType intervalSize=std::max<ScalarType>((ScalarType)UnifiedOrthoNode.size()/(ScalarType)SampleDensity,(ScalarType)1);
        ScalarType currPos=0;
        do
        {
            int currIndex=(int)floor(currPos+0.5);
            currIndex=std::min(currIndex,(int)UnifiedOrthoNode.size()-1);
            Nodes.push_back(UnifiedOrthoNode[currIndex]);
            currPos+=intervalSize;
        }while (currPos<UnifiedOrthoNode.size());
    }

    //    void GetOrthoGenerativeNodes(const std::vector<size_t> NodeI,
    //                                 std::vector<size_t> &OrthoNodes,
    //                                 size_t SampleDensity=20)
    //    {
    //        OrthoNodes.clear();
    //        std::vector<size_t> AllOrthoNodes,UnifiedOrthoNode;
    //        //AnisotropicQuery<MeshType>::getOrthogonalNodes(anigraph,NodeI,2,AllOrthoNodes);
    //        LoopFunctions<MeshType>::
    //        RetrieveCrossNodes(const AnisotropicGraph<MeshType> &anigraph,
    //                                           const std::vector<vcg::face::Pos<FaceType> > &FeaturePos,
    //                                           std::vector<size_t> &GuardNodes)
    //        //erase duplicates in position
    //        std::set<CoordType> NodePos;
    //        for (size_t i=0;i<AllOrthoNodes.size();i++)
    //        {
    //            CoordType currP=anigraph.NodePos(AllOrthoNodes[i]);
    //            if (NodePos.count(currP)>0)continue;
    //            NodePos.insert(currP);
    //            UnifiedOrthoNode.push_back(AllOrthoNodes[i]);
    //        }
    //        if (AllOrthoNodes.size()<SampleDensity)
    //        {
    //            OrthoNodes=AllOrthoNodes;
    //            return;
    //        }
    //        //then subsample
    //        ScalarType intervalSize=std::max<ScalarType>((ScalarType)UnifiedOrthoNode.size()/(ScalarType)SampleDensity,(ScalarType)1);
    //        ScalarType currPos=0;
    //        do
    //        {
    //            int currIndex=(int)floor(currPos+0.5);
    //            currIndex=std::min(currIndex,(int)UnifiedOrthoNode.size()-1);
    //            OrthoNodes.push_back(UnifiedOrthoNode[currIndex]);
    //            currPos+=intervalSize;
    //        }while (currPos<UnifiedOrthoNode.size());
    //        std::cout<<"Sampled "<<OrthoNodes.size()<<" out of "<<SampleDensity<<std::endl;
    //    }

    void GetOrthoGenerativeNodes(const std::vector<size_t> NodeI,
                                 std::vector<size_t> &OrthoNodes,
                                 size_t SampleDensity=20)
    {
        OrthoNodes.clear();
        AnisotropicQuery<MeshType>::getOrthogonalNodes(anigraph,NodeI,2,OrthoNodes);
        if (SampleDensity>0)
            SubsampleIntervals(OrthoNodes,SampleDensity);
    }

    void GetOrthoGenerativeNodes(const std::vector<vcg::face::Pos<FaceType> > &FeaturePos,
                                 std::vector<size_t> &OrthoNodes,
                                 size_t SampleDensity=20)
    {
        OrthoNodes.clear();
        //std::vector<size_t> AllOrthoNodes,UnifiedOrthoNode;
        //AnisotropicQuery<MeshType>::getOrthogonalNodes(anigraph,NodeI,2,AllOrthoNodes);
        LoopFunctions<MeshType>::RetrieveCrossNodes(anigraph,FeaturePos,OrthoNodes);
        SubsampleIntervals(OrthoNodes,SampleDensity);
        std::cout<<"Sampled "<<OrthoNodes.size()<<" out of "<<SampleDensity<<std::endl;
    }

    //    void GetOrthoGenerativeNodes(const std::vector<size_t> NodeI,
    //                                 std::vector<size_t> &OrthoNodes,
    //                                 ScalarType Lenght,
    //                                 size_t SampleDensity=5)
    //    {
    //        OrthoNodes.clear();

    //        std::vector<size_t> AllOrthoNodes;
    //        AnisotropicQuery<MeshType>::getOrthogonalNodes(anigraph,NodeI,2,AllOrthoNodes);

    //        //then subsample
    //        ScalarType intervalSize=(ScalarType)Lenght/(ScalarType)SampleDensity;
    //        CoordType currPos=anigraph.NodePos(AllOrthoNodes[0]);
    //        OrthoNodes.push_back(AllOrthoNodes[0]);
    //        for (size_t i=1;i<AllOrthoNodes.size();i++)
    //        {
    //            size_t currIndex=AllOrthoNodes[i];
    //            CoordType testPos=anigraph.NodePos(currIndex);
    //            if ((currPos-testPos).Norm()>=intervalSize)
    //            {
    //                currPos=testPos;
    //                OrthoNodes.push_back(currIndex);
    //            }
    //        }
    //    }

    //    void GetOrthoGenerativeNodes(const std::vector<size_t> NodeI,
    //                                 std::vector<size_t> &OrthoNodes,
    //                                 ScalarType Lenght,
    //                                 size_t SampleDensity=5)
    //    {
    //        OrthoNodes.clear();

    ////        std::vector<size_t> AllOrthoNodes;
    ////        AnisotropicQuery<MeshType>::getOrthogonalNodes(anigraph,NodeI,2,AllOrthoNodes);

    //        //then subsample
    //        ScalarType intervalSize=(ScalarType)Lenght/(ScalarType)SampleDensity;
    //        CoordType currPos=anigraph.NodePos(NodeI[0]);//AllOrthoNodes[0]);
    //        OrthoNodes.push_back(NodeI[0]);//AllOrthoNodes[0]);
    //        for (size_t i=1;i<NodeI.size();i++)
    //        {
    //            size_t currIndex=NodeI[i];//AllOrthoNodes[i];
    //            CoordType testPos=anigraph.NodePos(currIndex);
    //            if ((currPos-testPos).Norm()>=intervalSize)
    //            {
    //                currPos=testPos;
    //                OrthoNodes.push_back(currIndex);
    //            }
    //        }
    //    }

    void InitBarrierFromSharpFeatures(bool Reset=true)
    {
        if (Reset)
        {
            ActiveNodes.assign(ActiveNodes.size(),true);
        }

        //if there are already loops compute distances from them
        if (SharpFeatures.size()==0)return;


        for (size_t i=0;i<SharpFeatures.size();i++)
        {
            for (size_t j=0;j<SharpFeatures[i].TangentGuardNodes.size();j++)
            {
                int NodeId=SharpFeatures[i].TangentGuardNodes[j];
                ActiveNodes[NodeId]=false;
                int OppId=anigraph.OppositeNode(NodeId);
                ActiveNodes[OppId]=false;
            }
        }

        for (size_t i=0;i<SingNodes.size();i++)
        {
            int NodeId=SingNodes[i];
            ActiveNodes[NodeId]=false;
            int OppId=anigraph.OppositeNode(NodeId);
            ActiveNodes[OppId]=false;
        }

        SetBarrierNonActiveNodes();

    }

    void InitSharpFeatureIntersectingPos()
    {
        for (size_t i=0;i<SharpFeatures.size();i++)
            for (size_t j=0;j<SharpFeatures[i].GuardNodes.size();j++)
            {
                CoordType pos=anigraph.NodePos(SharpFeatures[i].GuardNodes[j]);
                SharpFeatures[i].IntersectingPos.insert(pos);
            }
    }

    std::vector<int> SharpIntersectionNeeds;
    std::vector<int> LoopIntersectionNeeds;

    std::vector<CoordType> Intersections;

    void GetIntersections(SharpFeature &Sf0,SharpFeature &Sf1,
                          std::vector<CoordType> &IntersectPos)
    {
        IntersectPos.clear();
        set_intersection(Sf0.IntersectingPos.begin(),Sf0.IntersectingPos.end(),
                         Sf1.IntersectingPos.begin(),Sf1.IntersectingPos.end(),
                         std::inserter(IntersectPos,IntersectPos.begin()));
    }

    void GetIntersections(SharpFeature &Sf,const CandidateLoop &Cl,
                          std::vector<CoordType> &IntersectPos)
    {
        IntersectPos.clear();
        std::set<CoordType> PosLoop;
        for (size_t i=0;i<Cl.LoopPath.nodes.size();i++)
        {
            size_t IndexN=Cl.LoopPath.nodes[i];
            PosLoop.insert(anigraph.NodePos(IndexN));
        }
        set_intersection(Sf.IntersectingPos.begin(),Sf.IntersectingPos.end(),
                         PosLoop.begin(),PosLoop.end(),
                         std::inserter(IntersectPos,IntersectPos.begin()));
    }

    void GetIntersections(const CandidateLoop &Cl0,const CandidateLoop &Cl1,
                          std::vector<CoordType> &IntersectPos)
    {
        IntersectPos.clear();
        ConflictFinder<MeshType> ConflFinder(anigraph);
        ConflFinder.FindCrossIntersections(Cl0.LoopPath,Cl1.LoopPath,true,true,IntersectPos);
    }

    void GetCandidateIntersections(const CandidateLoop &Cl,
                                   std::vector<int> &SharpInters,
                                   std::vector<CoordType> &SharpIntersPos,
                                   std::vector<int> &LoopInters,
                                   std::vector<CoordType> &LoopIntersPos)
    {
        for (int i=0;i<(int)SharpFeatures.size();i++)
        {
            std::vector<CoordType> IntersectPos;

            GetIntersections(SharpFeatures[i],Cl,IntersectPos);

            if (IntersectPos.size()==0)continue;

            //then add the info
            SharpInters.resize(SharpInters.size()+IntersectPos.size(),i);
            SharpIntersPos.insert(SharpIntersPos.end(),IntersectPos.begin(),IntersectPos.end());
        }

        for (int i=0;i<(int)ChoosenLoops.size();i++)
        {
            std::vector<CoordType> IntersectPos;

            GetIntersections(Cl,ChoosenLoops[i],IntersectPos);

            if (IntersectPos.size()==0)continue;

            LoopInters.resize(LoopInters.size()+IntersectPos.size(),i);
            LoopIntersPos.insert(LoopIntersPos.end(),IntersectPos.begin(),IntersectPos.end());
        }
    }

    //    void UpdateCrossIntersection()
    //    {
    //        std::cout<<"Check Sharp Sharp Intersections"<<std::endl;
    //        for (int i=0;i<(int)SharpFeatures.size()-1;i++)
    //        {
    //            if (!SharpFeatures[i].IsLoop)continue;
    //            for (int j=(i+1);j<(int)SharpFeatures.size();j++)
    //            {
    //                if (!SharpFeatures[j].IsLoop)continue;
    //                std::vector<CoordType> IntersectPos;
    //                GetIntersections(SharpFeatures[i],SharpFeatures[j],IntersectPos);
    //                if (IntersectPos.size()==0)continue;
    //                Intersections.insert(Intersections.end(),IntersectPos.begin(),IntersectPos.end());
    //                SharpIntersectionNeeds[i]-=IntersectPos.size();
    //                SharpIntersectionNeeds[j]-=IntersectPos.size();
    //                SharpIntersectionNeeds[i]=std::max(0,SharpIntersectionNeeds[i]);
    //                SharpIntersectionNeeds[j]=std::max(0,SharpIntersectionNeeds[j]);
    //            }
    //        }
    //        std::cout<<"done"<<std::endl;

    //        std::cout<<"Check Sharp Loop Intersections"<<std::endl;
    //        for (int i=0;i<(int)SharpFeatures.size();i++)
    //        {
    //            //if (!SharpFeatures[i].IsLoop)continue;
    //            for (int j=0;j<(int)ChoosenLoops.size();j++)
    //            {
    //                std::vector<CoordType> IntersectPos;

    //                GetIntersections(SharpFeatures[i],ChoosenLoops[j],IntersectPos);

    //                if (IntersectPos.size()==0)continue;
    //                Intersections.insert(Intersections.end(),IntersectPos.begin(),IntersectPos.end());
    //                SharpIntersectionNeeds[i]-=IntersectPos.size();

    //                if (SharpFeatures[i].IsConcave)
    //                    LoopIntersectionNeeds[j]-=IntersectPos.size();

    //                SharpIntersectionNeeds[i]=std::max(0,SharpIntersectionNeeds[i]);
    //                LoopIntersectionNeeds[j]=std::max(0,LoopIntersectionNeeds[j]);
    //            }
    //        }

    //        std::cout<<"done"<<std::endl;

    //        std::cout<<"Check Loops Loops Intersections"<<std::endl;
    //        for (int i=0;i<(int)ChoosenLoops.size()-1;i++)
    //            for (int j=(i+1);j<(int)ChoosenLoops.size();j++)
    //            {
    //                std::vector<CoordType> IntersectPos;
    //                GetIntersections(ChoosenLoops[i],ChoosenLoops[j],IntersectPos);
    //                if (IntersectPos.size()==0)continue;
    //                Intersections.insert(Intersections.end(),IntersectPos.begin(),IntersectPos.end());
    //                LoopIntersectionNeeds[i]-=IntersectPos.size();
    //                LoopIntersectionNeeds[j]-=IntersectPos.size();
    //                LoopIntersectionNeeds[i]=std::max(0,LoopIntersectionNeeds[i]);
    //                LoopIntersectionNeeds[j]=std::max(0,LoopIntersectionNeeds[j]);
    //            }
    //        std::cout<<"done"<<std::endl;
    //    }

    void AddIntersection(size_t Index0,size_t Index1,
                         bool Sharp0,bool Sharp1,
                         bool Cut0,bool Cut1)
    {
        if (Sharp0)
        {
            if (Cut1)
                Sharp_Cutting_Int[Index0]++;
            else
                Sharp_No_Cutting_Int[Index0]++;
        }

        if (Sharp1)
        {
            if (Cut0)
                Sharp_Cutting_Int[Index1]++;
            else
                Sharp_No_Cutting_Int[Index1]++;
        }

        if (!Sharp0)
        {
            if (Cut1)
                Loop_Cutting_Int[Index0]++;
            else
                Loop_No_Cutting_Int[Index0]++;
        }

        if (!Sharp1)
        {
            if (Cut0)
                Loop_Cutting_Int[Index1]++;
            else
                Loop_No_Cutting_Int[Index1]++;
        }
    }

    void UpdateIntersectionNeeds()
    {
        SharpIntersectionNeeds.clear();
        LoopIntersectionNeeds.clear();

        SharpIntersectionNeeds.resize(SharpFeatures.size(),Intersection_needs);
        LoopIntersectionNeeds.resize(ChoosenLoops.size(),Intersection_needs);

        assert(Sharp_Cutting_Int.size()==SharpFeatures.size());
        assert(Sharp_No_Cutting_Int.size()==SharpFeatures.size());
        assert(Loop_Cutting_Int.size()==ChoosenLoops.size());
        assert(Loop_No_Cutting_Int.size()==ChoosenLoops.size());

        for (size_t i=0;i<SharpFeatures.size();i++)
        {
            if (SharpFeatures[i].IsLoop)continue;
            SharpIntersectionNeeds[i]=0;
        }

        //remove the ones that intersect
        for (size_t i=0;i<SharpIntersectionNeeds.size();i++)
        {
            if (!SharpFeatures[i].IsLoop) continue;

            SharpIntersectionNeeds[i]-=Sharp_Cutting_Int[i];
            SharpIntersectionNeeds[i]-=Sharp_No_Cutting_Int[i];
        }

        //then set max to zero and check pairness of cutting ones
        for (size_t i=0;i<SharpIntersectionNeeds.size();i++)
        {
            if (!SharpFeatures[i].IsLoop) continue;

            SharpIntersectionNeeds[i]=std::max(SharpIntersectionNeeds[i],0);

            //            if ((SharpIntersectionNeeds[i]==0) && (Sharp_Cutting_Int[i]%2==1))
            //                SharpIntersectionNeeds[i]=1;
        }

        for (size_t i=0;i<LoopIntersectionNeeds.size();i++)
        {
            LoopIntersectionNeeds[i]-=Loop_Cutting_Int[i];
            LoopIntersectionNeeds[i]-=Loop_No_Cutting_Int[i];
        }

        //then set max to zero and check pairness of cutting ones
        for (size_t i=0;i<LoopIntersectionNeeds.size();i++)
        {
            LoopIntersectionNeeds[i]=std::max(LoopIntersectionNeeds[i],0);

            if ((LoopIntersectionNeeds[i]==0) && (Loop_Cutting_Int[i]==1))
                LoopIntersectionNeeds[i]=1;
        }
    }

    void UpdateLastLoopIntersection()
    {
        Loop_Cutting_Int.resize(ChoosenLoops.size(),0);
        Loop_No_Cutting_Int.resize(ChoosenLoops.size(),0);

        std::cout<<"Check Sharp Loop Intersections"<<std::endl;
        if (ChoosenLoops.size()==0)return;
        for (int i=0;i<(int)SharpFeatures.size();i++)
        {
            bool IsCutting0=SharpFeatures[i].IsConcave;
            IsCutting0&=SharpFeatures[i].IsLoop;
            bool IsCutting1=true;

            std::vector<CoordType> IntersectPos;

            GetIntersections(SharpFeatures[i],ChoosenLoops.back(),IntersectPos);

            if (IntersectPos.size()==0)continue;
            Intersections.insert(Intersections.end(),IntersectPos.begin(),IntersectPos.end());

            //std::cout<<"int "<<IntersectPos.size()<<std::endl;
            for (size_t k=0;k<IntersectPos.size();k++)
                AddIntersection(i,ChoosenLoops.size()-1,true,false,IsCutting0,IsCutting1);
        }
        if (ChoosenLoops.size()<2)return;
        for (int i=0;i<(int)ChoosenLoops.size()-1;i++)
        {
            bool IsCutting0=true;
            bool IsCutting1=true;
            std::vector<CoordType> IntersectPos;
            GetIntersections(ChoosenLoops[i],ChoosenLoops.back(),IntersectPos);
            if (IntersectPos.size()==0)continue;
            Intersections.insert(Intersections.end(),IntersectPos.begin(),IntersectPos.end());

            //std::cout<<"int "<<IntersectPos.size()<<std::endl;
            for (size_t k=0;k<IntersectPos.size();k++)
                AddIntersection(i,ChoosenLoops.size()-1,false,false,IsCutting0,IsCutting1);
        }
    }

    void UpdateCrossIntersection()
    {
        Sharp_Cutting_Int.clear();
        Sharp_Cutting_Int.resize(SharpFeatures.size(),0);

        Sharp_No_Cutting_Int.clear();
        Sharp_No_Cutting_Int.resize(SharpFeatures.size(),0);

        Loop_Cutting_Int.clear();
        Loop_Cutting_Int.resize(ChoosenLoops.size(),0);

        Loop_No_Cutting_Int.clear();
        Loop_No_Cutting_Int.resize(ChoosenLoops.size(),0);

        Intersections.clear();

        std::cout<<"Check Sharp Sharp Intersections"<<std::endl;
        for (int i=0;i<(int)SharpFeatures.size()-1;i++)
        {
            bool IsCutting0=SharpFeatures[i].IsConcave;
            for (int j=(i+1);j<(int)SharpFeatures.size();j++)
            {
                std::vector<CoordType> IntersectPos;
                GetIntersections(SharpFeatures[i],SharpFeatures[j],IntersectPos);
                if (IntersectPos.size()==0)continue;

                bool IsCutting1=SharpFeatures[j].IsConcave;

                //set the intersection point for debugging
                Intersections.insert(Intersections.end(),IntersectPos.begin(),IntersectPos.end());

                for (size_t k=0;k<IntersectPos.size();k++)
                    AddIntersection(i,j,true,true,IsCutting0,IsCutting1);
            }
        }
        std::cout<<"done"<<std::endl;

        std::cout<<"Check Sharp Loop Intersections"<<std::endl;
        if (ChoosenLoops.size()!=0)
        {
            for (int i=0;i<(int)SharpFeatures.size();i++)
            {
                bool IsCutting0=SharpFeatures[i].IsConcave;
                IsCutting0&=SharpFeatures[i].IsLoop;
                for (int j=0;j<(int)ChoosenLoops.size();j++)
                {
                    bool IsCutting1=true;

                    std::vector<CoordType> IntersectPos;

                    GetIntersections(SharpFeatures[i],ChoosenLoops[j],IntersectPos);

                    if (IntersectPos.size()==0)continue;
                    Intersections.insert(Intersections.end(),IntersectPos.begin(),IntersectPos.end());

                    //std::cout<<"int "<<IntersectPos.size()<<std::endl;
                    for (size_t k=0;k<IntersectPos.size();k++)
                        AddIntersection(i,j,true,false,IsCutting0,IsCutting1);
                }
            }
        }

        std::cout<<"done"<<std::endl;

        std::cout<<"Check Loops Loops Intersections"<<std::endl;
        if (ChoosenLoops.size()!=0)
        {
            for (int i=0;i<(int)ChoosenLoops.size()-1;i++)
                for (int j=(i+1);j<(int)ChoosenLoops.size();j++)
                {
                    bool IsCutting0=true;
                    bool IsCutting1=true;
                    std::vector<CoordType> IntersectPos;
                    GetIntersections(ChoosenLoops[i],ChoosenLoops[j],IntersectPos);
                    if (IntersectPos.size()==0)continue;
                    Intersections.insert(Intersections.end(),IntersectPos.begin(),IntersectPos.end());

                    //std::cout<<"int "<<IntersectPos.size()<<std::endl;
                    for (size_t k=0;k<IntersectPos.size();k++)
                        AddIntersection(i,j,false,false,IsCutting0,IsCutting1);
                }
        }
        std::cout<<"done"<<std::endl;
    }

    //    void InitIntersectionNeeds()
    //    {
    //        SharpIntersectionNeeds=std::vector<int> (SharpFeatures.size(),4);
    //        for (size_t i=0;i<SharpFeatures.size();i++)
    //        {
    //            if (SharpFeatures[i].IsLoop)continue;
    //            SharpIntersectionNeeds[i]=0;
    //        }
    //        LoopIntersectionNeeds=std::vector<int> (ChoosenLoops.size(),4);
    //    }


    void InitSharpFeaturesSimple()
    {

        //FeaturePos.clear();
        SharpFeatures.clear();

        //first find concave edges
        std::vector<std::vector<vcg::face::Pos<FaceType> > > ConcaveEdges;
        //anigraph.Mesh().GetSharpFeatures(ConcaveEdges,ConcaveEdge);
        sharp.GetSharpFeatures(ConcaveEdges,ConcaveEdge);

        //check if they are loops
        for (size_t i=0;i<ConcaveEdges.size();i++)
        {
            bool IsLoop=LoopFunctions<MeshType>::IsClosedLoop(ConcaveEdges[i]);
            SharpFeatures.push_back(SharpFeature(ConcaveEdges[i],IsLoop,true));
        }

        //then find convex edges
        std::vector<std::vector<vcg::face::Pos<FaceType> > > ConvexEdges;
        //anigraph.Mesh().GetSharpFeatures(ConvexEdges,ConvexEdge);
        sharp.GetSharpFeatures(ConvexEdges,ConvexEdge);

        //check if they are loops
        for (size_t i=0;i<ConvexEdges.size();i++)
        {
            bool IsLoop=LoopFunctions<MeshType>::IsClosedLoop(ConvexEdges[i]);
            SharpFeatures.push_back(SharpFeature(ConvexEdges[i],IsLoop,false));
        }

        //then initialize the guard nodes
        for (size_t i=0;i<SharpFeatures.size();i++)
        {
            //ScalarType LengthFeature=LoopFunctions<MeshType>::PosLenght(SharpFeatures[i].FeaturePos);
            LoopFunctions<MeshType>::RetrieveGuardNodes(anigraph,SharpFeatures[i].FeaturePos,SharpFeatures[i].TangentGuardNodes);

            if (SharpFeatures[i].IsLoop)
                GetOrthoGenerativeNodes(SharpFeatures[i].FeaturePos,SharpFeatures[i].OrthoGenerativeNodes,SParams.OrthoLoopsSample);

            //GetOrthoGenerativeNodes(SharpFeatures[i].GuardNodes,SharpFeatures[i].GenerativeNodes);

            std::vector<size_t> OrthoNodes;
            LoopFunctions<MeshType>::RetrieveCrossNodes(anigraph,SharpFeatures[i].FeaturePos,OrthoNodes);
            std::vector<CoordType> NodePositions;
            AnisotropicQuery<MeshType>::getNodesPositions(anigraph,OrthoNodes,NodePositions);
            SharpFeatures[i].IntersectingPos=std::set<CoordType>(NodePositions.begin(),NodePositions.end());
        }

        PrintFeaturesStats();

    }

    std::vector<PathType> ProblematicLoops;
    SharpFeatureSampler<MeshType> *SharpSampl;

    void InitSharpFeaturesTracing()
    {
        //        FeaturePos.clear();
        //        Concave.clear();

        std::vector<PathType> ConcaveLoops;
        std::vector<std::vector<vcg::face::Pos<FaceType> > > ClosedConcavePos, UnsampledConcavePos,ConvexPos;

        SharpSampl=new SharpFeatureSampler<MeshType>(anigraph,CTable,sharp);
        //SharpFeatureSampler<MeshType> SharpSampl(anigraph,CTable);

        //trace the loops
        //assert(0);
        //NEED FINISH


        SharpSampl->GetLoopsGreedyAdd(SParams.AvoidSelfIntersections,
                                      SParams.CloseConvexEndPoints,
                                      SParams.OnlyOrthogonalClose,
                                      ConcaveLoops,ProblematicLoops,
                                      ConvexPath,ClosedConcavePos,
                                      UnsampledConcavePos,ConvexPos,
                                      SParams.PenaltyDrift,SParams.MaxDev);


        for (size_t i=0;i<ConcaveLoops.size();i++)
        {

            //set basics data
            ChoosenLoops.push_back(CandidateLoop());
            ChoosenLoops.back().To_Update=false;
            ChoosenLoops.back().BunrdOut=false;
            ChoosenLoops.back().InitialNode=ConcaveLoops[i].nodes[0];
            ChoosenLoops.back().LoopPath=ConcaveLoops[i];
            ChoosenLoops.back().AvgDistance=0;

            //then set generative and tangent nodes
            LoopFunctions<MeshType>::GetTangentNodes(anigraph,ChoosenLoops.back().LoopPath,
                                                     ChoosenLoops.back().TangentGuardNodes);
            GetOrthoGenerativeNodes(ChoosenLoops.back().LoopPath.nodes,
                                    ChoosenLoops.back().OrthoGenerativeNodes);


        }

        //add concave closed from beginning
        for (size_t i=0;i<UnsampledConcavePos.size();i++)
        {
            bool IsLoop=LoopFunctions<MeshType>::IsClosedLoop(UnsampledConcavePos[i]);
            assert(!IsLoop);
            SharpFeatures.push_back(SharpFeature(UnsampledConcavePos[i],IsLoop,true));
        }


        for (size_t i=0;i<ClosedConcavePos.size();i++)
        {
            bool IsLoop=LoopFunctions<MeshType>::IsClosedLoop(ClosedConcavePos[i]);
            assert(IsLoop);
            SharpFeatures.push_back(SharpFeature(ClosedConcavePos[i],IsLoop,true));
        }

        for (size_t i=0;i<ConvexPos.size();i++)
        {
            bool IsLoop=LoopFunctions<MeshType>::IsClosedLoop(ConvexPos[i]);
            SharpFeatures.push_back(SharpFeature(ConvexPos[i],IsLoop,false));
        }

        //then initialize the guard nodes
        for (size_t i=0;i<SharpFeatures.size();i++)
        {
            //ScalarType LengthFeature=LoopFunctions<MeshType>::PosLenght(SharpFeatures[i].FeaturePos);
            LoopFunctions<MeshType>::RetrieveGuardNodes(anigraph,SharpFeatures[i].FeaturePos,SharpFeatures[i].TangentGuardNodes);

            if (SharpFeatures[i].IsLoop)
                GetOrthoGenerativeNodes(SharpFeatures[i].FeaturePos,SharpFeatures[i].OrthoGenerativeNodes);

            //GetOrthoGenerativeNodes(SharpFeatures[i].GuardNodes,SharpFeatures[i].GenerativeNodes);

            std::vector<size_t> OrthoNodes;
            LoopFunctions<MeshType>::RetrieveCrossNodes(anigraph,SharpFeatures[i].FeaturePos,OrthoNodes);
            std::vector<CoordType> NodePositions;
            AnisotropicQuery<MeshType>::getNodesPositions(anigraph,OrthoNodes,NodePositions);
            SharpFeatures[i].IntersectingPos=std::set<CoordType>(NodePositions.begin(),NodePositions.end());
        }
        //std::cout<<"E"<<std::endl;
    }

    //    void UpdateDistByLoop(CandidateLoop &CLoop)
    //    {
    //        std::vector<PathType> TestPaths;
    //        TestPaths.push_back(CLoop.LoopPath);
    //        size_t t0=clock();

    //        std::cout<<"*** Update distances considering new added loop ***"<<std::endl;

    //        //        PathType testP;
    //        //        for (size_t i=0;i<SingNodes.size();i++)
    //        //        {
    //        //            if (!ActiveNodes[SingNodes[i]])continue;
    //        //            testP.nodes.push_back(SingNodes[i]);
    //        //        }
    //        //        TestPaths.push_back(testP);

    //        //std::vector<size_t> singNodes,sinkNodes;
    //        //        anigraph.GetSingularityNodes(singNodes,sinkNodes);
    //        //        PathType testP;
    //        //        std::cout<<"Sing Nodes "<<singNodes.size()<<std::endl;
    //        //        for (size_t i=0;i<singNodes.size();i++)
    //        //            testP.nodes.push_back(singNodes[i]);
    //        //        TestPaths.push_back(testP);
    //        //NodeDist.clear();

    //        AnisotropicQuery<MeshType>::ComputeLoopDistances(anigraph,TestPaths,45,0,NodeDist);
    //        size_t t1=clock();
    //        TimePriorityUpdate+=t1-t0;
    //        std::cout<<"*** Done ***"<<std::endl;
    //    }

    void GetCandidateGuardNodes(CandidateLoop &CLoop,std::vector<int> &GuardNodes)
    {
        GuardNodes.clear();
        for (size_t i=0;i<CLoop.TangentGuardNodes.size();i++)
        {
            GuardNodes.insert(GuardNodes.end(),
                              CLoop.TangentGuardNodes[i].begin(),
                              CLoop.TangentGuardNodes[i].end());
        }
    }

    void GetChoosenLoopsGuardNodes(std::vector<int> &GuardNodes)
    {
        GuardNodes.clear();
        for (size_t i=0;i<ChoosenLoops.size();i++)
        {
            std::vector<int> GuardNodesTemp;
            GetCandidateGuardNodes(ChoosenLoops[i],GuardNodesTemp);
            GuardNodes.insert(GuardNodes.end(),GuardNodesTemp.begin(),GuardNodesTemp.end());
        }
    }

    void UpdateDistByLoop(CandidateLoop &CLoop)
    {
        //        std::vector<PathType> TestPaths;
        //        TestPaths.push_back(CLoop.LoopPath);
        size_t t0=clock();

        std::cout<<"*** Update distances considering new added loop ***"<<std::endl;
        std::vector<int> GuardNodes;
        GetCandidateGuardNodes(CLoop,GuardNodes);
        AnisotropicQuery<MeshType>::ComputeDistanceFromNodes(anigraph,GuardNodes,45,0,NodeDist);
        size_t t1=clock();
        TimePriorityUpdate+=t1-t0;
        std::cout<<"*** Done ***"<<std::endl;
    }

    void UpdateDistChoosenLoops()
    {
        //        std::vector<PathType> ChoosenPaths;
        //        GetChoosenPaths(ChoosenPaths);

        size_t t0=clock();
        std::cout<<"*** Update distances from choosen loops ***"<<std::endl;
        std::vector<int> GuardNodes;
        GetChoosenLoopsGuardNodes(GuardNodes);
        AnisotropicQuery<MeshType>::ComputeDistanceFromNodes(anigraph,GuardNodes,45,0,NodeDist);
        //AnisotropicQuery<MeshType>::ComputeDistanceFromLoops(anigraph,ChoosenPaths,45,0,NodeDist);
        size_t t1=clock();
        TimePriorityUpdate+=t1-t0;
        std::cout<<"*** Done ***"<<std::endl;
    }

    void UpdateDistSharpFeatures()
    {
        //        std::vector<PathType> ChoosenPaths;
        //        GetChoosenPaths(ChoosenPaths);

        size_t t0=clock();

        NodeDist.clear();




        std::cout<<"*** Update distances from features loops ***"<<std::endl;
        std::vector<size_t> SharpNodes;
        for (size_t i=0;i<SharpFeatures.size();i++)
        {
            //std::cout<<"*** De ***"<<std::endl;
            SharpNodes.insert(SharpNodes.end(),SharpFeatures[i].TangentGuardNodes.begin(),SharpFeatures[i].TangentGuardNodes.end());
        }
        //        //PathType testP;
        //        for (size_t i=0;i<SingNodes.size();i++)
        //        {
        //            if (!ActiveNodes[SingNodes[i]])continue;
        //            //testP.nodes.push_back(SingNodes[i]);
        //            SharpNodes.push_back(SingNodes[i]);
        //        }
        //SharpNodes.push_back(testP);

        AnisotropicQuery<MeshType>::ComputeNodesDistances(anigraph,SharpNodes,45,0,NodeDist);

        size_t t1=clock();
        TimePriorityUpdate+=t1-t0;
        std::cout<<"*** Done ***"<<std::endl;
    }

    void InitSharpFeatures()
    {
        if (!SParams.FindConcaveLoops)
            InitSharpFeaturesSimple();
        else
            InitSharpFeaturesTracing();

    }

    void AddCandidatesFromSharpFeature(size_t IndexFeature)
    {
        for (size_t i=0;i<SharpFeatures[IndexFeature].OrthoGenerativeNodes.size();i++)
        {
            int CurrNode=SharpFeatures[IndexFeature].OrthoGenerativeNodes[i];

            CandidateLoop CLoop(CurrNode);
            Candidates.push_back(CLoop);
        }
    }

    void AddCandidatesFromSharpFeatures()
    {
        std::cout<<"*** Adding Candidates from sharp features ***"<<std::endl;
        for (size_t i=0;i<SharpFeatures.size();i++)
        {
            if (!SharpFeatures[i].IsLoop)continue;
            AddCandidatesFromSharpFeature(i);
        }
        std::cout<<"*** Done ***"<<std::endl;
    }

    std::vector<CoordType> GraphOriginalPos;

    void SaveOriginalPos()
    {
        GraphOriginalPos.clear();
        for (size_t i=0;i<anigraph.Graph.size();i++)
            GraphOriginalPos.push_back(anigraph.Graph[i]->pos);
    }

    void InitStats()
    {
        TimeFirstLoops=0;
        TimeFirstGuess=0;
        TimeLoopUpdate=0;
        TimePriorityUpdate=0;
        InitialLoopNum=0;
        UpdatedLoopNum=0;
    }

    void Reset()
    {
        Candidates.clear();
        ChoosenLoops.clear();
        NodeDist.clear();
        ActiveNodes=std::vector<bool>(anigraph.NumNodes(),true);
        NodeDist=std::vector<ScalarType>(anigraph.NumNodes(),std::numeric_limits<ScalarType>::max());
        Intersections.clear();

        Sharp_Cutting_Int.clear();
        Sharp_No_Cutting_Int.clear();
        Loop_Cutting_Int.clear();
        Loop_No_Cutting_Int.clear();
    }


    bool UpdatedDist;

    //    void InitVertQBySingDistance()
    //    {
    //        std::vector<VertexType*> SingVert;

    //        for (size_t i=0;i<anigraph.Mesh().vert.size();i++)
    //        {
    //            bool IsSing=vcg::tri::CrossField<MeshType>::IsSingular(anigraph.Mesh(),anigraph.Mesh().vert[i]);
    //            if (!IsSing)continue;
    //            SingVert.push_back(&anigraph.Mesh().vert[i]);
    //        }

    //        vcg::tri::Geodesic<MeshType>::Compute(anigraph.Mesh(),SingVert);
    //        ScalarType MaxQ=0;
    //        for (size_t i=0;i<anigraph.Mesh().vert.size();i++)
    //            MaxQ=std::max(MaxQ,anigraph.Mesh().vert[i].Q());

    //        //clamp
    //        //MaxQ*=0.05;
    //        for (size_t i=0;i<anigraph.Mesh().vert.size();i++)
    //        {
    //            anigraph.Mesh().vert[i].Q()=std::min(anigraph.Mesh().vert[i].Q(),MaxQ);
    //            anigraph.Mesh().vert[i].Q()/=MaxQ;
    //            //anigraph.Mesh().vert[i].Q()=0.0001+(1-anigraph.Mesh().vert[i].Q());
    //            //if (anigraph.Mesh().vert[i].Q()<)
    //            //anigraph.Mesh().vert[i].Q()=pow(anigraph.Mesh().vert[i].Q(),2);
    //        }
    //        //anigraph.WeightLenghtByVertexQ();
    //    }

    ConflictTable<MeshType> CTable;

public:

    void RestoreOriginalPos()
    {
        for (size_t i=0;i<anigraph.Graph.size();i++)
            anigraph.Graph[i]->pos=GraphOriginalPos[i];
    }

    struct SampleParams
    {
        ScalarType MaxDev;
        ScalarType PenaltyDrift;
        //ScalarType SharpDegree;
        //ScalarType FlatDegree;
        int AdditionalLoops;
        bool FindConcaveLoops;
        bool AvoidSelfIntersections;
        bool ResampleAdditionalLoops;
        bool UpdateLoops;
        bool CloseConvexEndPoints;
        bool OnlyOrthogonalClose;
        bool FixInterCross;
        int OrthoLoopsSample;
        int OversampleRate;

        SampleParams()
        {
            MaxDev=45;
            PenaltyDrift=100;
            OversampleRate=20;
            //SharpDegree=60;
            FindConcaveLoops=true;
            AvoidSelfIntersections=true;
            ResampleAdditionalLoops=true;
            UpdateLoops=true;
            CloseConvexEndPoints=true;
            OnlyOrthogonalClose=true;
            AdditionalLoops=100;
            FixInterCross=true;
            OrthoLoopsSample=20;
        }
    };

    SampleParams SParams;

    void UpdateLoopsStep()
    {
        //update the loops
        UpdateCandidates();
    }

    //    void InitSingNodes()
    //    {
    //        SingNodes.clear();

    //        std::set<size_t> CrossFaces;
    //        for (size_t i=0;i<anigraph.Mesh().face.size();i++)
    //        {
    //            bool HasSing=false;
    //            for (size_t j=0;j<3;j++)
    //            {
    //                VertexType *v=anigraph.Mesh().face[i].V(j);
    //                HasSing|=vcg::tri::CrossField<MeshType>::IsSingular(anigraph.Mesh(),*v);
    //            }
    //            if (HasSing)
    //                CrossFaces.insert(i);
    //        }

    //        for (size_t i=0;i<anigraph.NumNodes();i++)
    //        {
    //            if (!anigraph.IsEdgeNode(i))continue;
    //            size_t IndexF,M4Dir;
    //            anigraph.GetFaceDir(i,IndexF,M4Dir);
    //            if (CrossFaces.count(IndexF)==0)continue;
    //            SingNodes.push_back(i);
    //            size_t IndexOpp=anigraph.OppositeNode(i);
    //            SingNodes.push_back(IndexOpp);
    //        }
    //        SingNodesSet=std::set<size_t>(SingNodes.begin(),SingNodes.end());
    //    }

    void InitSampling(SampleParams &_SParams)
    {

        SParams=_SParams;
        InitStats();

        Reset();

        //InitSingNodes()
        //        if (GraphOriginalPos.size()==anigraph.NumNodes())
        //            RestoreOriginalPos();



        InitSharpFeatures();

        //TODO if there is not enough then add samples poisson
        //InitIntersectionNeeds();

        //add the candidates from the sharp features
        AddCandidatesFromSharpFeatures();
        UpdatedDist=false;


        anigraph.SetAllValid();

        if (SharpFeatures.size()>0)
        {
            UpdateDistSharpFeatures();
            UpdatedDist=true;
        }
        if (ChoosenLoops.size()>0)
        {
            UpdateDistChoosenLoops();
            UpdatedDist=true;
        }

        //set barrier from sharp features
        InitBarrierFromSharpFeatures();


        for (size_t i=0;i<ChoosenLoops.size();i++)
        {
            std::vector<int> ConflNodes;
            CTable.GetConflictNodes(ChoosenLoops[i].LoopPath,ConflNodes);
            for (size_t j=0;j<ConflNodes.size();j++)
            {
                assert(ConflNodes[j]>=0);
                assert(ConflNodes[j]<(int)ActiveNodes.size());
                ActiveNodes[ConflNodes[j]]=false;
            }
        }

        for (size_t i=0;i<ConvexPath.size();i++)
        {
            std::vector<int> ConflNodes;
            CTable.GetConflictNodes(ConvexPath[i],ConflNodes);
            for (size_t j=0;j<ConflNodes.size();j++)
            {
                assert(ConflNodes[j]>=0);
                assert(ConflNodes[j]<(int)ActiveNodes.size());
                ActiveNodes[ConflNodes[j]]=false;
            }
        }

        SetBarrierNonActiveNodes();
        for (size_t i=0;i<ChoosenLoops.size();i++)
            CTable.InvalidateTangantEdges(anigraph,ChoosenLoops[i].LoopPath);

        for (size_t i=0;i<ConvexPath.size();i++)
            CTable.InvalidateTangantEdges(anigraph,ConvexPath[i]);


        UpdateLoopsStep();
    }

    void AddLoopsFrom(CandidateLoop &CLoop,bool SubSample=true)
    {

        if (SubSample)
            GetOrthoGenerativeNodes(CLoop.LoopPath.nodes,CLoop.OrthoGenerativeNodes,SParams.OrthoLoopsSample);
        else
            GetOrthoGenerativeNodes(CLoop.LoopPath.nodes,CLoop.OrthoGenerativeNodes,0);
        //GetOrthoGenerativeNodes(CLoop,CLoop.GenerativeNodes);
        std::cout<<"Found "<<CLoop.OrthoGenerativeNodes.size()<<" generative nodes"<<std::endl;
        for (size_t i=0;i<CLoop.OrthoGenerativeNodes.size();i++)
        {
            int IndexNode=CLoop.OrthoGenerativeNodes[i];
            if (!ActiveNodes[IndexNode])continue;
            CandidateLoop newCLoop(IndexNode);
            Candidates.push_back(newCLoop);
        }
    }

    void PrintStack()
    {
        for (size_t i=0;i<Candidates.size();i++)
            std::cout<<"Solve Sharp "<<Candidates[i].SolveSharp<<" Solve Loop "<<Candidates[i].SolveLoop<<" Dist "<<Candidates[i].AvgDistance<<std::endl;
    }

    std::vector<size_t> Sharp_Cutting_Int;
    std::vector<size_t> Sharp_No_Cutting_Int;
    std::vector<size_t> Loop_Cutting_Int;
    std::vector<size_t> Loop_No_Cutting_Int;

    //    void SampleStep()
    //    {
    //        if (Candidates.size()==0)return;
    //        if ((Candidates.back().SolveSharp==0)&&(Candidates.back().SolveLoop==0))return;
    //        PrintStack();

    //        //then add one candidate
    //        assert(Candidates.back().IsValid());
    //        ChoosenLoops.push_back(Candidates.back());
    //        //add new needs for intersections

    //        Candidates.pop_back();

    //        //then add new candidates
    //        //        GetOrthoGenerativeNodes(SharpFeatures[i].GuardNodes,SharpFeatures[i].GenerativeNodes);
    //        //        ChoosenLoops.back()
    //        //        //Update barrier of sharp features
    //        //        InitBarrierFromSharpFeatures();

    //        //then update the distances considering the new added
    //        UpdateDistByLoop(ChoosenLoops.back());

    //        //update barrier considering last one
    //        std::vector<PathType> TestLoops;
    //        TestLoops.push_back(ChoosenLoops.back().LoopPath);

    //        //add the new need
    //        LoopIntersectionNeeds.push_back(ChoosenLoops.back().OpenNum);

    //        std::vector<bool> IsLoop(1,true);

    //        std::vector<size_t> ToDisableNodes;
    //        LoopFunctions<MeshType>::BarrierPaths(anigraph,TestLoops,IsLoop,&ToDisableNodes);
    //        for (size_t j=0;j<ToDisableNodes.size();j++)
    //        {
    //            //ChoosenLoops.back().GuardNodes.push_back(ToDisableNodes[j]);
    //            ActiveNodes[ToDisableNodes[j]]=false;
    //        }

    //        //get the ones need to be recomputed
    //        for (size_t i=0;i<Candidates.size();i++)
    //            SetUpdateFlag(Candidates[i]);

    //        //set new seeds
    //        if (SParams.ResampleAdditionalLoops)
    //            AddLoopsFrom(ChoosenLoops.back());

    //        //then recompute new loops
    //        UpdateLoopsStep();

    //        //add the intersections
    //        CandidateLoop InsertedLoop=ChoosenLoops.back();
    //        Intersections.insert(Intersections.end(),InsertedLoop.SharpIntersPos.begin(),
    //                             InsertedLoop.SharpIntersPos.end());
    //        Intersections.insert(Intersections.end(),InsertedLoop.LoopIntersPos.begin(),
    //                             InsertedLoop.LoopIntersPos.end());

    //        //update current needs
    //        for (size_t i=0;i<InsertedLoop.SharpInters.size();i++)
    //        {
    //            int IndexShF=InsertedLoop.SharpInters[i];
    //            SharpIntersectionNeeds[IndexShF]--;
    //            SharpIntersectionNeeds[IndexShF]=std::max(0,SharpIntersectionNeeds[IndexShF]);
    //        }

    //        for (size_t i=0;i<InsertedLoop.LoopInters.size();i++)
    //        {
    //            int IndexL=InsertedLoop.LoopInters[i];
    //            LoopIntersectionNeeds[IndexL]--;
    //            LoopIntersectionNeeds[IndexL]=std::max(0,LoopIntersectionNeeds[IndexL]);
    //        }

    //        //then compute the priority
    //        UpdateCandidatePriorityByLoopDistances();

    //        //update priority by intersections
    //        UpdateCandidatePriorityByIntersections();

    //        //then sort the candidates
    //        std::sort(Candidates.begin(),Candidates.end());

    //        //PrintStack();
    //        //        //add the new need
    //        //        LoopIntersectionNeeds.push_back(InsertedLoop.OpenNum);

    //        //UpdateCrossIntersection();
    //    }
    std::vector<int> ProblematicNodes;

    void UpdateBarrierFromLoop(CandidateLoop &CLoop)
    {
        std::vector<int> ConflNodes;
        CTable.GetConflictNodes(CLoop.LoopPath,ConflNodes);
        for (size_t j=0;j<ConflNodes.size();j++)
        {
            assert(ConflNodes[j]>=0);
            assert(ConflNodes[j]<(int)ActiveNodes.size());
            ActiveNodes[ConflNodes[j]]=false;
        }
        SetBarrierNonActiveNodes();
        for (size_t i=0;i<ChoosenLoops.size();i++)
            CTable.InvalidateTangantEdges(anigraph,ChoosenLoops[i].LoopPath);

        for (size_t i=0;i<ConvexPath.size();i++)
            CTable.InvalidateTangantEdges(anigraph,ConvexPath[i],false);
    }

    bool SampleStep(bool essential=true)
    {
        //PrintStack();

        if (Candidates.size()==0)return false;
        //if (Candidates.back().SolveSharp==0)return false;

        if ((essential)&&(Candidates.back().SolveSharp==0)&&(Candidates.back().SolveLoop==0))return false;

        //then add one candidate
        assert(Candidates.back().IsValid());
        ChoosenLoops.push_back(Candidates.back());
        ChoosenLoops.back().IsEssential=essential;
        //add new needs for intersections

        Candidates.pop_back();

        //then update the distances considering the new added
        UpdateDistByLoop(ChoosenLoops.back());

        //update self distances
        //UpdateCandidatesSelfDistances();

        //update barrier considering last one
        //        std::vector<PathType> TestLoops;
        //        TestLoops.push_back(ChoosenLoops.back().LoopPath);

        //disable nodes
        UpdateBarrierFromLoop(ChoosenLoops.back());

        //        std::vector<int> ConflNodes;
        //        CTable.GetConflictNodes(ChoosenLoops.back().LoopPath,ConflNodes);
        //        for (size_t j=0;j<ConflNodes.size();j++)
        //        {
        //            assert(ConflNodes[j]>=0);
        //            assert(ConflNodes[j]<ActiveNodes.size());
        //            ActiveNodes[ConflNodes[j]]=false;
        //        }
        //        SetBarrierNonActiveNodes();
        //        for (size_t i=0;i<ChoosenLoops.size();i++)
        //            CTable.InvalidateTangantEdges(anigraph,ChoosenLoops[i].LoopPath);

        //        for (size_t i=0;i<ConvexPath.size();i++)
        //            CTable.InvalidateTangantEdges(anigraph,ConvexPath[i],false);


        //        std::vector<bool> IsLoop(1,true);
        //        std::vector<size_t> ToDisableNodes;
        //        LoopFunctions<MeshType>::BarrierPaths(anigraph,TestLoops,IsLoop,&ToDisableNodes);
        //        for (size_t j=0;j<ToDisableNodes.size();j++)
        //        {
        //            //ChoosenLoops.back().GuardNodes.push_back(ToDisableNodes[j]);
        //            ActiveNodes[ToDisableNodes[j]]=false;
        //        }

        //get the ones need to be recomputed
        for (size_t i=0;i<Candidates.size();i++)
            SetUpdateFlag(Candidates[i]);

        //set new seeds
        if (SParams.ResampleAdditionalLoops)
            (ChoosenLoops.back());

        //then recompute new loops
        UpdateLoopsStep();


        //update intersections
        std::cout<<"Updating Cross Intersections"<<std::endl;
        //UpdateCrossIntersection();
        UpdateLastLoopIntersection();
        std::cout<<"done"<<std::endl;
        std::cout<<"Updating Intersection needs"<<std::endl;
        UpdateIntersectionNeeds();
        std::cout<<"done"<<std::endl;

        //then compute the priority
        std::cout<<"Updating Distances"<<std::endl;
        UpdateCandidatePriorityByLoopDistances();
        std::cout<<"done"<<std::endl;

        //update priority by intersections
        std::cout<<"Updating Intersections"<<std::endl;
        UpdateCandidatePriorityByIntersections();
        std::cout<<"done"<<std::endl;

        //then sort the candidates
        std::sort(Candidates.begin(),Candidates.end());


        //UpdateDistByLoop(ChoosenLoops.back());

        //        ProblematicNodes.clear();
        //        std::vector<int> TestNodes;
        //        GetCandidateGuardNodes(Candidates.back(),TestNodes);
        //        //GetCandidateGuardNodes(ChoosenLoops.back(),TestNodes);
        //        for (size_t i=0;i<TestNodes.size();i++)
        //        {
        //            int IndexN=TestNodes[i];
        //            if(!ActiveNodes[IndexN])
        //            {
        //            ProblematicNodes.push_back(IndexN);
        //            std::cout<<"Error "<<IndexN<<std::endl;
        //            }
        //            //assert(anigraph.IsValid(IndexN));
        //        }

        //for (sizeCandidates.back()
        //PrintStack();
        //        //add the new need
        //        LoopIntersectionNeeds.push_back(InsertedLoop.OpenNum);

        //UpdateCrossIntersection();
        return true;
    }


    //    bool FixInterCrossStep(size_t InsertPos)
    //    {
    //        InterCrossValidator<MeshType> ICross(anigraph);
    //        std::vector<PathType> ChoosenPaths;
    //        GetChoosenPaths(ChoosenPaths);
    //        std::vector<typename InterCrossValidator<MeshType>::UnsolvedInfo> CrossUnsolved;

    //        std::cout<<"1- Gathering Unsolved ones"<<std::endl;
    //        ICross.GetCrossUnsolved(ChoosenPaths,CrossUnsolved);
    //        size_t UnsolvedBefore=CrossUnsolved.size();
    //        if (UnsolvedBefore==0)return false;

    //        std::cout<<"2- Gathering nodes to trace from"<<std::endl;
    //        Candidates.clear();
    //        std::vector<size_t> To_Sample;
    //        std::set<CoordType> StartPoint;
    //        for (size_t i=0;i<CrossUnsolved.size();i++)
    //        {
    //            size_t IndexP=CrossUnsolved[i].PathIndex;
    //            size_t IndexN0=CrossUnsolved[i].Start;
    //            size_t IndexN1=CrossUnsolved[i].End;
    //            size_t SizeLoop=ChoosenPaths[IndexP].nodes.size();
    //            for (size_t i=IndexN0;i!=IndexN1;i=(i+1)%SizeLoop)
    //            {
    //                size_t IndexN=ChoosenPaths[IndexP].nodes[i];
    //                CoordType Pos=anigraph.NodePos(IndexN);
    //                StartPoint.insert(Pos);
    //            }
    //            To_Sample.push_back(IndexP);
    //        }

    //        std::sort(To_Sample.begin(),To_Sample.end());
    //        To_Sample.erase( std::unique( To_Sample.begin(), To_Sample.end() ), To_Sample.end() );

    //        for (size_t i=0;i<To_Sample.size();i++)
    //        {
    //            AddLoopsFrom(ChoosenLoops[To_Sample[i]],false);
    //        }
    //        //filter the ones not in the map of pos
    //        std::vector<CandidateLoop> CandidatesSwap;
    //        std::cout<<"Candidate before fileter "<<Candidates.size()<<std::endl;

    //        for (size_t i=0;i<Candidates.size();i++)
    //        {
    //            size_t IndexN=Candidates[i].InitialNode;
    //            CoordType Pos=anigraph.NodePos(IndexN);
    //            if (StartPoint.count(Pos)==0)continue;
    //            CandidatesSwap.push_back(Candidates[i]);
    //        }
    //        Candidates=CandidatesSwap;
    //        std::cout<<"Candidate TO trace "<<Candidates.size()<<std::endl;
    //        //update the  loops
    //        UpdateLoopsStep();
    //        std::cout<<"Candidates traced "<<Candidates.size()<<std::endl;

    //        //then update priorities
    //        UpdateCandidatePriorityByLoopDistances();

    //        std::cout<<"Checking Inter Cross"<<std::endl;
    //        CandidatesSwap.clear();
    //        for (size_t i=0;i<Candidates.size();i++)
    //        {
    //            //add one further loop
    //            ChoosenPaths.push_back(Candidates[i].LoopPath);
    //            ICross.GetCrossUnsolved(ChoosenPaths,CrossUnsolved,false);
    //            ChoosenPaths.pop_back();
    //            size_t UnsolvedAfter=CrossUnsolved.size();
    //            int NumSolved=UnsolvedBefore-UnsolvedAfter;
    //            if (NumSolved<=0)continue;
    //            //            ChoosenLoops.push_back(Candidates[i]);
    //            Candidates[i].SolveInterCross=NumSolved;
    //            CandidatesSwap.push_back(Candidates[i]);
    //            //            std::cout<<"Can Solve "<<CandidatesSwap.back().SolveInterCross<<std::endl;
    //        }
    //        Candidates=CandidatesSwap;
    //        if (Candidates.size()==0)return false;

    //        std::sort(Candidates.begin(),Candidates.end());

    //        //if (Candidates.back().SolveInterCross==0)return false;

    //        ChoosenLoops.insert(ChoosenLoops.begin()+InsertPos,Candidates.back());
    //        std::cout<<"Solved Loops "<<Candidates.back().SolveInterCross<<std::endl;
    //        UpdateBarrierFromLoop(ChoosenLoops.back());
    //        return true;
    //    }

    bool TestAddingLoop(std::vector<std::vector<vcg::face::Pos<FaceType> > > &Features,
                        const InterCrossValidator<MeshType> &ICross,
                        std::vector<PathType> &ChoosenPaths,
                        const size_t &UnsolvedBefore,
                        const size_t &TangentNode,
                        CandidateLoop &NewPath)
    {
        assert(TangentNode<anigraph.NumNodes());

        //get the orthogonal one (need to change the name of the function)
        std::vector<size_t> OrthoN;
        anigraph.TangentialNodes(TangentNode,OrthoN,2);
        if (OrthoN.size()==0)return false;
        for (size_t i=0;i<OrthoN.size();i++)
            if (anigraph.IsValid(OrthoN[i]))
            {
                Candidates.clear();
                NewPath=CandidateLoop(OrthoN[i]);
                Candidates.push_back(NewPath);

                //then trace add one as candidate
                UpdateLoopsStep();

                //if not traced
                if (Candidates.size()==0)continue;

                //get back the traced
                NewPath=Candidates.back();

                //check if solved the issue
                ChoosenPaths.push_back(Candidates.back().LoopPath);
                std::vector<typename InterCrossValidator<MeshType>::UnsolvedInfo> CrossUnsolved;
                ICross.GetCrossUnsolved(Features,ChoosenPaths,CrossUnsolved,false);
                size_t UnsolvedAfter=CrossUnsolved.size();
                int NumSolved=UnsolvedBefore-UnsolvedAfter;
                if (NumSolved>0)
                {
                    ChoosenPaths.pop_back();
                    return true;
                }
                else
                    ChoosenPaths.pop_back();
            }
        return false;
    }

    bool TestAddingLoop(std::vector<std::vector<vcg::face::Pos<FaceType> > > &Features,
                        const InterCrossValidator<MeshType> &ICross,
                        std::vector<PathType> &ChoosenPaths,
                        const size_t &UnsolvedBefore,
                        const std::vector<size_t> &TangentNodes,
                        CandidateLoop &NewPath)
    {
        int Interval=(TangentNodes.size()/SParams.OrthoLoopsSample);
        Interval=std::max(Interval,1);

        for (size_t i=0;i<TangentNodes.size();i+=Interval)
        {
            bool Solved=TestAddingLoop(Features,ICross,ChoosenPaths,UnsolvedBefore,TangentNodes[i],NewPath);
            if (!Solved)continue;
            return true;
        }
        return false;
    }

    bool FixInterCrossStep(std::vector<std::vector<vcg::face::Pos<FaceType> > > &Features,
                           size_t InsertPos)
    {
        InterCrossValidator<MeshType> ICross(anigraph);
        std::vector<PathType> ChoosenPaths;
        GetChoosenPaths(ChoosenPaths);
        std::vector<typename InterCrossValidator<MeshType>::UnsolvedInfo> CrossUnsolved;

        std::cout<<"1- Gathering Unsolved ones"<<std::endl;
        //        std::vector<std::vector<vcg::face::Pos<FaceType> > > Features;
        //        for (size_t i=0;i<SharpFeatures.size();i++)
        //            Features.push_back(SharpFeatures[i].FeaturePos);

        ICross.GetCrossUnsolved(Features,ChoosenPaths,CrossUnsolved,false);
        size_t UnsolvedBefore=CrossUnsolved.size();
        if (UnsolvedBefore==0)return false;
//        for (size_t i=0;i<ChoosenPaths.size();i++)
//            std::cout<<"Test "<<ChoosenPaths[i].nodes.size()<<std::endl;

        std::cout<<"2- start trying to solve "<<std::endl;
        for (size_t i=0;i<CrossUnsolved.size();i++)
        {
            //go over each loop
            size_t IndexN0=CrossUnsolved[i].Start;
            size_t IndexN1=CrossUnsolved[i].End;

            std::vector<size_t> TangentNodes;
            size_t IndexP=CrossUnsolved[i].PathIndex;
            assert(IndexP<ChoosenPaths.size());
            size_t SizeLoop=ChoosenPaths[IndexP].nodes.size();
//            std::cout<<"Index "<<IndexP<<std::endl;
//            std::cout<<"Start "<<IndexN0<<std::endl;
//            std::cout<<"End "<<IndexN1<<std::endl;
//            std::cout<<"Size "<<SizeLoop<<std::endl;
            assert(IndexN0<SizeLoop);
            assert(IndexN1<SizeLoop);

            for (size_t i=IndexN0;i!=IndexN1;i=(i+1)%SizeLoop)
            {
                size_t IndexN=ChoosenPaths[IndexP].nodes[i];
                TangentNodes.push_back(IndexN);
            }
            //std::cout<<"- 3 "<<std::endl;
            CandidateLoop NewPath;
            //std::cout<<"- test adding "<<std::endl;
            if (!TestAddingLoop(Features,ICross,ChoosenPaths,UnsolvedBefore,TangentNodes,NewPath))continue;

            ChoosenLoops.insert(ChoosenLoops.begin()+InsertPos,NewPath);
            std::cout<<"ADDED LOOP"<<std::endl;
            //std::cout<<"Solved Loops "<<Candidates.back().SolveInterCross<<std::endl;
            UpdateBarrierFromLoop(ChoosenLoops[InsertPos]);
            return true;
        }

        return false;
    }

    bool CleanRemainingConflicts()
    {
        std::vector<PathType> ChoosenPaths;
        GetChoosenPaths(ChoosenPaths);
        ConflictFinder<MeshType> ConflFinder(anigraph);
        std::vector<ConflictInfo> Conflicts;
        ConflFinder.FindConflicts(ChoosenPaths,Conflicts,true,true);
        if (Conflicts.size()==0)return false;

        std::cout<<"***WARNING "<<Conflicts.size()<<" CONFLICTS ***"<<std::endl;
        std::cout<<"...Cleaning..."<<std::endl;
        std::vector<std::pair<size_t,size_t> > ConflNumber;
        for (size_t i=0;i<ChoosenPaths.size();i++)
            ConflNumber.push_back(std::pair<size_t,size_t>(0,i));

        for (size_t i=0;i<Conflicts.size();i++)
        {
            std::cout<<"*PATHS: ";

            for (size_t j=0;j<Conflicts[i].Paths.size();j++)
            {
                size_t Index=Conflicts[i].Paths[j];
                ConflNumber[Index].first++;
                std::cout<<" "<<Index;
            }
            assert(Conflicts[i].Paths.size()>=2);
            std::cout<<std::endl;
        }
        std::sort(ConflNumber.begin(),ConflNumber.end());
        //at least one conflict stored
        assert(ConflNumber.back().first>0);
        size_t ToRemoveIndex=ConflNumber.back().second;
        DeleteLoop(ToRemoveIndex);
        return true;
    }

    void SampleLoops(SampleParams &_SParams)
    {
        //LoopFunctions<MeshType>::CheckOppositeValid(anigraph);

        std::cout<<"Initializing Intersections Table"<<std::endl;
        CTable.Init(anigraph);
        std::cout<<"done"<<std::endl;

        //return;
        //InitVertQBySingDistance();

        InitSampling(_SParams);


        UpdateCrossIntersection();

        UpdateIntersectionNeeds();

        //then compute the priority
        UpdateCandidatePriorityByLoopDistances();

        //update priority by intersections
        UpdateCandidatePriorityByIntersections();

        //then sort the candidates
        std::sort(Candidates.begin(),Candidates.end());

        AddingEssential=true;

        std::cout<<"Sampling Essential Loops"<<std::endl;
        while (SampleStep()){}

        size_t InsertPos=ChoosenLoops.size();

        AddingEssential=false;

        //if (SParams.AdditionalLoops==0)return;

        //disable sharp features faces
        //LoopFunctions<MeshType>::DisableFaceSharp(anigraph,anigraph.Mesh().FeatureSeq);
        LoopFunctions<MeshType>::DisableFaceSharp(anigraph,sharp.FeatureSeq);
        for (size_t i=0;i<anigraph.NumNodes();i++)
        {
            if (anigraph.IsValid(i))continue;
            ActiveNodes[i]=false;
        }

        std::cout<<"Sampling Additional Loops"<<std::endl;
        int MinSamples=std::max(50,SParams.AdditionalLoops*SParams.OversampleRate);
        SampleCandidates(MinSamples);
        //then recompute new loops
        UpdateLoopsStep();
        //if ((ChoosenLoops.size()>0)&&())
        if (UpdatedDist)
            UpdateCandidatePriorityByLoopDistances();
        else
        {
            std::cout<<"updating priority by self distance"<<std::endl;
            UpdatePriorityBySelfDist();
        }

        //then sort the candidates
        std::sort(Candidates.begin(),Candidates.end());

        size_t Loops0=ChoosenLoops.size();
        bool terminated=(SParams.AdditionalLoops==0);
        if (terminated)return;
        while (SampleStep(false)&&(!terminated))
        {
            size_t Loops1=ChoosenLoops.size();
            terminated=((int)(Loops1-Loops0)>=SParams.AdditionalLoops);
        }
        std::cout<<"Added "<<ChoosenLoops.size()-Loops0 <<" Additional Loops"<<std::endl;

//        for (size_t i=0;i<ChoosenLoops.size();i++)
//            std::cout<<"Test 0"<<ChoosenLoops[i].LoopPath.nodes.size()<<std::endl;


        std::vector<std::vector<vcg::face::Pos<FaceType> > > Features;
        for (size_t i=0;i<SharpFeatures.size();i++)
            Features.push_back(SharpFeatures[i].FeaturePos);
        if (SParams.FixInterCross)
        {
            //            FixInterCrossStep(InsertPos);
            //            FixInterCrossStep(InsertPos);
            //            FixInterCrossStep(InsertPos);
            std::cout<<"*** FIXING UNSOLVED PATHS ***"<<std::endl;
            while (FixInterCrossStep(Features,InsertPos)){};

            UpdateCrossIntersection();
            UpdateIntersectionNeeds();
        }

        InterCrossValidator<MeshType> ICross(anigraph);
        std::vector<PathType> ChoosenPaths;
        GetChoosenPaths(ChoosenPaths);
        std::vector<typename InterCrossValidator<MeshType>::UnsolvedInfo> CrossUnsolved;
        ICross.GetCrossUnsolved(Features,ChoosenPaths,CrossUnsolved,false);
        if (CrossUnsolved.size()>0)
        {
            std::cout<<"*** UNSOLVED PATHS ***"<<std::endl;
            for (size_t i=0;i<CrossUnsolved.size();i++)
            {
                size_t IndexP=CrossUnsolved[i].PathIndex;
                std::cout<<"Path Index "<<IndexP<<std::endl;

                LoopIntersectionNeeds[IndexP]=1;
            }
        }
        while(CleanRemainingConflicts()){};
    }

    void GetLoopSmoothPos(const PathType &p,std::vector<CoordType> &SPos)
    {
        SPos.clear();
        int sizeL=p.nodes.size();
        SPos.resize(p.nodes.size(),CoordType(0,0,0));
        for (int j=0;j<sizeL;j++)
        {
            size_t Nj=p.nodes[j];
            size_t N0=p.nodes[(j+sizeL-1)%sizeL];
            size_t N1=p.nodes[(j+1)%sizeL];
            CoordType Pos0=anigraph.NodePos(N0);
            CoordType Pos1=anigraph.NodePos(N1);
            SPos[j]=(Pos0+Pos1)/2;

            size_t FaceId,EdgeId,M4Dir;
            anigraph.GetFaceEdgeDir(Nj,FaceId,EdgeId,M4Dir);

            CoordType PosE1=anigraph.Mesh().face[FaceId].P0(EdgeId);
            CoordType PosE2=anigraph.Mesh().face[FaceId].P1(EdgeId);
            vcg::Segment3<ScalarType> s(PosE1,PosE2);
            ScalarType dist;
            vcg::SegmentPointDistance(s,SPos[j],SPos[j],dist);
        }
    }

    void SnapLoop()
    {
        std::set<size_t> SnappedNodes;
        //first fidn possible problems
        if (SharpSampl!=NULL)
            SharpSampl->SnapLoops(SnappedNodes);

        if (SnappedNodes.size()==0)return;

        //save all possible nodes that might have problems
        //move the nodes that can intersect after the corner in order to void the
        //intersection to collapse on a point
        std::map<std::pair<CoordType,CoordType>,std::vector<size_t> > CheckNodes;
        std::map<size_t,CoordType> SmoothedPos;
        for (size_t i=0;i<ChoosenLoops.size();i++)
        {
            int sizeL=ChoosenLoops[i].LoopPath.nodes.size();
            std::vector<CoordType> SPos;
            GetLoopSmoothPos(ChoosenLoops[i].LoopPath,SPos);

            for (int j=0;j<sizeL;j++)
            {
                size_t N0=ChoosenLoops[i].LoopPath.nodes[j];
                size_t N1=ChoosenLoops[i].LoopPath.nodes[(j+1)%sizeL];
                bool InN0=SnappedNodes.count(N0)>0;
                bool InN1=SnappedNodes.count(N1)>0;
                if (InN0 && InN1)continue;
                if ((!InN0) && (!InN1))continue;
                //check when only one of the two nodes is
                //snapped and move the other to the smooth position
                CoordType Pos0=anigraph.NodePos(N0);
                CoordType Pos1=anigraph.NodePos(N1);
                std::pair<CoordType,CoordType> key(std::min(Pos0,Pos1),std::max(Pos0,Pos1));
                if (InN0)
                {
                    CheckNodes[key].push_back(N1);
                    SmoothedPos[N1]=SPos[(j+1)%sizeL];
                }
                if (InN1)
                {
                    CheckNodes[key].push_back(N0);
                    SmoothedPos[N0]=SPos[j];
                }
            }
        }
        //then find possible conflicts
        typename std::map<std::pair<CoordType,CoordType>,std::vector<size_t> >::iterator IteM;
        for (IteM=CheckNodes.begin();IteM!=CheckNodes.end();IteM++)
        {
            if ((*IteM).second.size()>1)
            {
                std::cout<<"*** WARNING: ONE NODE CLOSE TO INTERSECTION THAT NEED TO BE DISPLACED ***"<<std::endl;
                assert((*IteM).second.size()<3);
                //move the second (it is the same)
                size_t ToMove0=(*IteM).second[0];
                size_t ToMove1=(*IteM).second[1];
                assert(SmoothedPos.count(ToMove0)>0);
                assert(SmoothedPos.count(ToMove1)>0);
                CoordType POld0=anigraph.NodePos(ToMove0);
                CoordType POld1=anigraph.NodePos(ToMove1);
                CoordType PNew0=SmoothedPos[ToMove0];
                CoordType PNew1=SmoothedPos[ToMove1];
                anigraph.Graph[ToMove0]->pos=POld0*0.5+PNew0*0.5;
                anigraph.Graph[ToMove1]->pos=POld1*0.5+PNew1*0.5;
                //               size_t FaceId,EdgeId,M4Dir;
                //               anigraph.GetFaceEdgeDir(ToMove,FaceId,EdgeId,M4Dir);
                //               CoordType Pos0=anigraph.NodePos(ToMove);
                //               CoordType Pos1=anigraph.Mesh().face[FaceId].P0(EdgeId);
                //               CoordType Pos2=anigraph.Mesh().face[FaceId].P1(EdgeId);
                //               ScalarType Dist1=(Pos0-Pos1).Norm();
                //               ScalarType Dist2=(Pos0-Pos2).Norm();
                //               CoordType TargetPos=Pos0*0.7+Pos1*0.3;
                //               if (Dist2<Dist1)
                //                   TargetPos=Pos0*0.7+Pos2*0.3;
                //               //then finally move!
                //               anigraph.Graph[ToMove]->pos=TargetPos;
            }
        }
    }

#ifndef NO_TRACING_OPENGL
    //**** DRAWING FUNCTIONS ****
    void GlDrawDisabledLinks(vcg::Color4b currC=vcg::Color4b(255,0,0,255))
    {
        glPushAttrib(GL_ALL_ATTRIB_BITS);
        glDisable(GL_LIGHTING);
        glDepthRange(0,0.99);


        glPointSize(10);
        vcg::glColor(vcg::Color4b(255,0,255,255));
        glBegin(GL_POINTS);

        for (size_t i=0;i<ProblematicNodes.size();i++)
        {
            CoordType Pos=anigraph.NodePos(ProblematicNodes[i]);
            vcg::glVertex(Pos);
            //anigraph.GlDrawNodeDebugInfo(currN,45,2,0.001);
        }

        glEnd();


        glLineWidth(2);

        vcg::glColor(currC);

        glBegin(GL_LINES);

        for (size_t Index0=0;Index0<anigraph.Graph.size();Index0++)
        {
            vcg::Point3f pos0;
            pos0.Import<ScalarType>(anigraph.Graph[Index0]->pos);
            for (typename AnisotropicGraph<MeshType>::NeighborsIterator nIt = anigraph.Graph[Index0]->neighbors.begin();
                 nIt != anigraph.Graph[Index0]->neighbors.end(); ++nIt)
            {
                if (nIt->Valid)continue;

                //get the neighbors node
                size_t Index1=nIt->nodeID;

                vcg::Point3f pos1;
                pos1.Import<ScalarType>(anigraph.Graph[Index1]->pos);

                vcg::glVertex(pos0);
                vcg::glVertex(pos1);
            }

        }
        glEnd();
        glPopAttrib();
    }

    void DrawGenerativeNodes()
    {
        glPushAttrib(GL_ALL_ATTRIB_BITS);
        glDisable(GL_LIGHTING);
        glDepthRange(0,0.999);
        glPointSize(10);
        vcg::glColor(vcg::Color4b(0,255,0,255));
        glBegin(GL_POINTS);

        for (size_t i=0;i<SharpFeatures.size();i++)
            for (size_t j=0;j<SharpFeatures[i].GenerativeNodes.size();j++)
            {
                size_t currN=SharpFeatures[i].GenerativeNodes[j];
                CoordType Pos=anigraph.NodePos(currN);
                vcg::glVertex(Pos);
                //anigraph.GlDrawNodeDebugInfo(currN,45,2,0.001);
            }
        for (size_t i=0;i<ChoosenLoops.size();i++)
            for (size_t j=0;j<ChoosenLoops[i].GenerativeNodes.size();j++)
            {
                size_t currN=ChoosenLoops[i].GenerativeNodes[j];
                CoordType Pos=anigraph.NodePos(currN);
                vcg::glVertex(Pos);
                //anigraph.GlDrawNodeDebugInfo(currN,45,2,0.001);
            }
        glEnd();
        glPopAttrib();
    }

    void DrawIntersections()
    {
        glPushAttrib(GL_ALL_ATTRIB_BITS);
        glDisable(GL_LIGHTING);
        glDepthRange(0,0.999);
        glPointSize(30);
        vcg::glColor(vcg::Color4b(255,0,0,255));
        glBegin(GL_POINTS);

        for (size_t i=0;i<Intersections.size();i++)
        {
            CoordType Pos=Intersections[i];
            vcg::glVertex(Pos);
            //anigraph.GlDrawNodeDebugInfo(currN,45,2,0.001);
        }
        glEnd();
        glPopAttrib();
    }

    void DrawFeatures()
    {

        glPushAttrib(GL_ALL_ATTRIB_BITS);
        glDisable(GL_LIGHTING);
        glDepthRange(0,0.999);
        glLineWidth(10);
        for (size_t i=0;i<SharpFeatures.size();i++)
        {

            vcg::Color4b currC=vcg::Color4b::Blue;
            if (SharpIntersectionNeeds[i]>0)
                currC=vcg::Color4b::Red;

            vcg::glColor(currC);
            glBegin(GL_LINES);
            for (size_t j=0;j<SharpFeatures[i].FeaturePos.size();j++)
            {
                vcg::face::Pos<FaceType> CurrPos=SharpFeatures[i].FeaturePos[j];
                CoordType P0=CurrPos.V()->P();
                CurrPos.FlipV();
                CoordType P1=CurrPos.V()->P();
                vcg::glVertex(P0);
                vcg::glVertex(P1);
            }
            glEnd();
        }

        glPopAttrib();
    }

    void DrawPath(const PathType &currP,vcg::Color4b path_col,
                  bool loop)
    {

        vcg::glColor(path_col);
        glPushAttrib(GL_ALL_ATTRIB_BITS);
        glDisable(GL_LIGHTING);

        glDepthRange(0,0.99998);
        glLineWidth(12);
        glBegin(GL_LINE_STRIP);
        for (size_t i=0;i<currP.nodes.size();++i)
        {
            vcg::glColor(path_col);
            vcg::glVertex(anigraph.NodePos(currP.nodes[i]));
        }
        if (loop)
            vcg::glVertex(anigraph.NodePos(currP.nodes[0]));
        glEnd();

        glDepthRange(0,0.99999);
        glLineWidth(16);
        glBegin(GL_LINE_STRIP);
        for (size_t i=0;i<currP.nodes.size();++i)
        {
            vcg::glColor(vcg::Color4b(0,0,0,255));
            vcg::glVertex(anigraph.NodePos(currP.nodes[i]));
        }
        if (loop)
            vcg::glVertex(anigraph.NodePos(currP.nodes[0]));
        glEnd();


        //        glPointSize(10);
        //        glBegin(GL_POINTS);
        //        for (size_t i=0;i<currP.nodes.size();++i)
        //        {
        //            vcg::glVertex(anigraph.NodePos(currP.nodes[i]));
        //        }
        //        glEnd();
        glPopAttrib();
    }

    //draw all loops
    void DrawLoops()
    {

        if (ChoosenLoops.size()==0)return;
        for (size_t i=0;i<ChoosenLoops.size();i++)
        {
            vcg::Color4b currC=vcg::Color4b::Blue;
            if (LoopIntersectionNeeds[i]>0)
                currC=vcg::Color4b::Red;

            if (selected==(int)i)
                currC=vcg::Color4b::Magenta;
            if (ChoosenLoops[i].LType==GuardLoop)continue;

            //vcg::Color4b currC=vcg::Color4b::Scatter(20,i%20,0.8);

            DrawPath(ChoosenLoops[i].LoopPath,currC,true);
            //anigraph.GlDrawPath(ChoosenLoops[i].LoopPath,currC,30,true);
        }

        for (size_t i=0;i<ConvexPath.size();i++)
        {
            vcg::Color4b currC=vcg::Color4b::Green;

            //anigraph.GlDrawPath(ConvexPath[i],currC,20,false);
            DrawPath(ChoosenLoops[i].LoopPath,currC,true);
        }

        for (size_t i=0;i<ProblematicLoops.size();i++)
        {
            vcg::Color4b currC=vcg::Color4b::Yellow;

            anigraph.GlDrawPath(ProblematicLoops[i],currC,4,true);
        }

        if (selected_loop==-1)return;

        //std::cout<<"AAA"<<std::endl;
        glPushAttrib(GL_ALL_ATTRIB_BITS);
        glDisable(GL_LIGHTING);
        glDepthRange(0,0.999);
        glPointSize(30);
        vcg::glColor(vcg::Color4b(0,0,0,255));
        int IndexNode=ChoosenLoops[selected_loop].LoopPath.nodes[selected_node];
        CoordType Pos=anigraph.NodePos(IndexNode);
        glBegin(GL_POINTS);
        vcg::glVertex(Pos);
        glEnd();
        glPopAttrib();
    }

    void DrawSingNodes()
    {
        for (size_t i=0;i<SingNodes.size();i++)
            anigraph.GlDrawNodeDebugInfo(SingNodes[i],45,1,0.001,true,false);

    }

    //draw candidates
    void DrawCandidates()
    {
        if (Candidates.size()==0)return;
        for (size_t i=0;i<Candidates.size();i++)
        {
            //vcg::Color4b currC=vcg::Color4b::ColorRamp(0,Candidates.size(),i);
            vcg::Color4b currC;
            currC.SetColorRampParula(0,Candidates.size(),i);
            //if (ChoosenLoops[i].LType==GuardLoop)continue;

            anigraph.GlDrawPath(Candidates[i].LoopPath,currC,2,true);
        }
    }
#endif

    void SelectClosestLoop(CoordType pos)
    {
        std::vector<PathType> testPath;
        for (size_t i=0;i<ChoosenLoops.size();i++)
            testPath.push_back(ChoosenLoops[i].LoopPath);

        for (size_t i=0;i<ConvexPath.size();i++)
            testPath.push_back(ConvexPath[i]);

        selected=LoopFunctions<MeshType>::GetClosestLoop(anigraph,testPath,pos);
        if (selected>=(int)ChoosenLoops.size())
        {
            ConvexPath.erase(ConvexPath.begin()+(selected-ChoosenLoops.size()));
            selected=-1;
        }
    }

    void SelectClosestLoopNode(CoordType pos)
    {
        std::vector<PathType> testPath;
        for (size_t i=0;i<ChoosenLoops.size();i++)
            testPath.push_back(ChoosenLoops[i].LoopPath);

        LoopFunctions<MeshType>::GetClosestLoopNode(anigraph,testPath,pos,selected_loop,selected_node);
        std::cout<<"Loop "<<selected_loop<<" Index "<<selected_node<<std::endl;
    }

    int GetClosestNodeGraph(CoordType pos)
    {
        ScalarType minD=std::numeric_limits<ScalarType>::max();
        int ClosestID=-1;
        for (size_t i=0;i<anigraph.NumNodes();i++)
        {
            CoordType test_pos=anigraph.NodePos(i);
            if ((test_pos-pos).Norm()>minD)continue;
            minD=(test_pos-pos).Norm();
            ClosestID=i;
        }
        return(ClosestID);
    }

    void MoveSelectedToNodeGraph(CoordType pos)
    {
        int ClosNode=GetClosestNodeGraph(pos);
        CoordType posN=anigraph.NodePos(ClosNode);
        ScalarType CurrD=(pos-posN).Norm();

        //find also the closest vert
        for (size_t i=0;i<anigraph.Mesh().vert.size();i++)
        {
            ScalarType TestD=(pos-anigraph.Mesh().vert[i].P()).Norm();
            if (TestD>CurrD)continue;
            posN=anigraph.Mesh().vert[i].P();
        }
        int IndexNode=ChoosenLoops[selected_loop].LoopPath.nodes[selected_node];
        anigraph.Graph[IndexNode]->pos=posN;//anigraph.NodePos(ClosNode);
    }


    void DeleteLoop(int selected)
    {
        ChoosenLoops.erase(ChoosenLoops.begin()+selected);
        LoopIntersectionNeeds.erase(LoopIntersectionNeeds.begin()+selected);

    }

    //    void UpdateCurvatureLoops(int RingSize)
    //    {
    //        vcg::tri::UP
    //    }
};


#endif //LOOP_OPTIMIZER
