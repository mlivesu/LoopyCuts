#ifndef SEPARATRIX_OPTIMIZER
#define SEPARATRIX_OPTIMIZER

//#include <anisotropic_geodesic.h>
#include <tracing_field/graph/anisotropic_graph.h>
#include <tracing_field/conflict_finder.h>
#include "scip_solver.h"
#include <tracing_field/stats_collector.h>
#include <edge_mesh_type.h>
#include "tracing_field/anisotropic_geodesic.h"

enum SepWeightType{SWLenght,SWUniform,SWMix};
enum SolveMode{SMConstraintSol,SMIterativeConstr,SMItaretiveGlobal};

/**
 * @brief The class that performs parametrization by considering the connections within the graph
 */
template < class MeshType >
class SeparatrixOptimizer {
    typedef typename MeshType::CoordType    CoordType;
    typedef typename MeshType::VertexType   VertexType;
    typedef typename MeshType::FaceType     FaceType;
    typedef typename MeshType::ScalarType   ScalarType;

public:
    /**
     * @brief AnisotropicGraph Constructor.
     * @param _anigraph
     */
    SeparatrixOptimizer(AnisotropicGraph<MeshType> &_anigraph,
                        StatsCollector<ScalarType> &_Stats) : anigraph(_anigraph),Stats(_Stats)
    {
        PercentileCutOff=false;
        MaxIterations=10;
    }

    /**
     * @brief ~AnisotropicGraph Destructor.
     */
    ~SeparatrixOptimizer() {}

private:


    AnisotropicGraph<MeshType> &anigraph;               //the mesh on which the separatrix graph is calculated
    StatsCollector<ScalarType> &Stats;

    typedef typename AnisotropicGraph<MeshType>::Path PathType;
    typedef typename ConflictFinder<MeshType>::ConflictInfo ConflictInfo;

    std::vector<size_t> InitSingNodes;                     //the vector of initial singularity nodes

public:
    //vector of found solutions
    std::vector<PathType> TotalSolution;
    bool PercentileCutOff;
    ScalarType cutoffL;
    MyEMesh CandidateEMesh;

private:

    //this class maintain all info per solve cycle
    std::vector<size_t> CurrSingNodes;                      //the vector of unsatisfied nodex
    std::vector<std::pair<size_t,size_t> > CurrSeparatrices;//the set of separatrices as two index connecting singularity nodes
    std::vector<PathType> CurrPaths;                        //the vector of paths
    std::vector<PathType> AllCandidates;                    //the vector of all candidates
    std::vector<ConflictInfo> CurrConflicts;                //face and paths that generate the conflict

    std::set<std::pair<int,int> > AddedPath;
    std::map<size_t,size_t> SingToIndex;                    //For each Singularity Node the variable index

    std::vector<size_t> NumEvaluated;                       //the number of evaluated separatrix for each singularity
    std::vector<size_t> TargetEvaluated;                    //the number of separatrix that have to be evaluated in total

    ScalarType minL,maxL;

    std::vector<size_t> UnsolvedNodes;                      //the vector of unsolved Nodes

    void Initialize()
    {
        //get the singularities and sinks informations
        std::vector<std::vector<size_t > > SingNodes,SinkNodes;
        std::vector<size_t> SingVertex,ExpValence;
        anigraph.GetSingularitiesInfo(SingNodes,SinkNodes,SingVertex,ExpValence);


        //then initialize the vector of initial singularities
        for (size_t i=0;i<SingNodes.size();i++)
            for (size_t j=0;j<SingNodes[i].size();j++)
                InitSingNodes.push_back(SingNodes[i][j]);

        //the vector of current singularities is equal to all singularities
        CurrSingNodes= std::vector<size_t>(InitSingNodes.begin(),InitSingNodes.end());
    }

    void AddPath(const PathType &CurrPath)
    {
        //get the two singularities
        //which are the first and the last of the path
        size_t SingNode0=CurrPath.nodes.front();
        size_t SinkNode1=CurrPath.nodes.back();

        //check the first must be a singulatity
        assert(anigraph.IsSingularity(SingNode0));
        //check the last one must be a sink
        assert(anigraph.IsSink(SinkNode1));

        //only insert in one sense the index of the nodes MUST be different
        assert(SingNode0!=SinkNode1);
        size_t SingNode1=anigraph.OppositeSink(SinkNode1);
        assert(anigraph.IsSingularity(SingNode1));

        //if start and end in the same singularity node then stop
        if (SingNode0==SingNode1)return;

        std::pair<int,int> key(std::min(SingNode0,SingNode1),std::max(SingNode0,SingNode1));
        if (AddedPath.count(key)>0)return;

        //check if self intersect
        ConflictFinder<MeshType> CFinder(anigraph);
        if (CFinder.SelfConflict(CurrPath,2,false))return;

        AddedPath.insert(key);

        //        //get the index of the singularities (can be also equal)
        //        int SingV0=anigraph.GetVertexSing(NodeSep0);
        //        int SingV1=anigraph.GetVertexSing(NodeSep1);

        //then get the indexes in the solver
        assert(SingToIndex.count(SingNode0)>0);
        assert(SingToIndex.count(SingNode1)>0);

        //get the two indices of the singularities
        int IndexVar0=SingToIndex[SingNode0];
        int IndexVar1=SingToIndex[SingNode1];

        //then add the path
        CurrPaths.push_back(CurrPath);
        AllCandidates.push_back(CurrPath);

        //and the separatrix for the system
        CurrSeparatrices.push_back(std::pair<size_t,size_t>(IndexVar0,IndexVar1));

        //update max and min lenght for visualization purposes
        ScalarType PathL=CurrPath.length;
        if (PathL<minL)minL=PathL;
        if (PathL>maxL)maxL=PathL;
    }

    void AddPaths(const std::vector<PathType> &TempPaths)
    {
        for (size_t i=0;i<TempPaths.size();i++)
        {
            //check the path
            assert(TempPaths[i].nodes.size()>0);
            AddPath(TempPaths[i]);
        }
    }

    bool FindSeparatrixStep(size_t maxConn,
                            ScalarType MaxDev,
                            ScalarType PenaltyDrift,
                            bool progressivePaths)
    {
        int t0 = clock();
        std::cout << "Computing Paths to singularities..." << std::endl;

        if (!progressivePaths)
        {
            CurrPaths.clear();
            CurrSeparatrices.clear();
            //reset the set of added to avoid multiple insertions
            AddedPath.clear();
        }

        //get the current step for visualize the progress
        size_t division=10;
        size_t step=floor((ScalarType)CurrSingNodes.size()/((ScalarType)division)+0.5)+1;

        //initialize the mapping
        //reset the structure that keep singularity
        SingToIndex.clear();
        for (size_t i=0;i<CurrSingNodes.size();i++)
        {
            size_t NodeSing=CurrSingNodes[i];
            SingToIndex[NodeSing]=i;
        }

        //then for each singularity find the possible connections to sinks
        for (size_t i=0;i<CurrSingNodes.size();i++)
        {
            //associate the index to the node
            size_t NodeSing=CurrSingNodes[i];
            //SingToIndex[NodeSing]=i;

            //safety check that is a singularity
            assert(anigraph.IsSingularity(NodeSing));

            //get all path emanating from the singularity
            std::vector<PathType> TempPaths;
            ScalarType MaxL=std::numeric_limits<ScalarType>::max();
            if (PercentileCutOff)MaxL=cutoffL;
            if (!progressivePaths)
                //anigraph.FindPatternToSingularities(NodeSing,MaxDev,PenaltyDrift,TempPaths,maxConn,false,MaxL);
                AnisotropicQuery<MeshType>::FindPatternToSingularities(anigraph,NodeSing,MaxDev,PenaltyDrift,TempPaths,maxConn,false,MaxL);
            else
                if (NumEvaluated[i]<TargetEvaluated[i])
                {
                    //anigraph.FindPatternToSingularities(NodeSing,MaxDev,PenaltyDrift,TempPaths,TargetEvaluated[i],false,MaxL);
                    AnisotropicQuery<MeshType>::FindPatternToSingularities(anigraph,NodeSing,MaxDev,PenaltyDrift,TempPaths,TargetEvaluated[i],false,MaxL);
                    NumEvaluated[i]=TargetEvaluated[i];
                }
            //print the progress
            ScalarType ratio=((ScalarType)i)/((ScalarType)CurrSingNodes.size());
            int progress=floor(ratio*100+0.5);
            if ((i%step)==0) std::cout << progress <<" %" << std::endl;

            //safety check paths should be in the right order
            if (TempPaths.size()==0)continue;

            //sefety check
            for (size_t i=0;i<TempPaths.size()-1;i++)
                assert(TempPaths[i].length<=TempPaths[i+1].length);

            AddPaths(TempPaths);
        }

        return (CurrPaths.size()>0);
    }

    void SetSolutionAsBarrier(bool checkReal)
    {
        ConflictFinder<MeshType> ConflFinder(anigraph);
        ConflFinder.SetGeodesicBarrier(TotalSolution,checkReal);
    }

    void FindObjectiveFunct(bool ConstrainSolution,
                            const SepWeightType &SepWType,
                            std::vector<ScalarType> &ObjFunct)
    {
        ObjFunct.clear();
        assert(CurrSeparatrices.size()==CurrPaths.size());
        for (size_t i=0;i<CurrPaths.size();i++)
        {
            ScalarType currW=1;
            if ((SepWType==SWLenght)||
                    (SepWType==SWMix))
                currW=CurrPaths[i].length;

            if (ConstrainSolution)
                ObjFunct.push_back(currW);
            else
            {
                if (SepWType==SWMix)
                    ObjFunct.push_back(1+1.0/currW);
                else
                    ObjFunct.push_back(1.0/currW);
            }
        }
    }

    void FindSingularityConstraints(std::vector<std::vector<ScalarType> > &CoeffsV,
                                    std::vector<ScalarType> &LhsV,
                                    std::vector<ScalarType> &RhsV,
                                    bool ConstrainSolution)
    {
        CoeffsV.clear();
        LhsV.clear();
        RhsV.clear();

        //resize: 1 constraint per Singularity Node and 1 variable for each path
        CoeffsV.resize(CurrSingNodes.size(),std::vector<ScalarType>(CurrPaths.size(),0));

        //set the intervals of the constraint
        if (!ConstrainSolution)
            LhsV.resize(CurrSingNodes.size(),0);
        else
            LhsV.resize(CurrSingNodes.size(),1);

        RhsV.resize(CurrSingNodes.size(),1);

        //safety check
        assert(CurrPaths.size()==CurrPaths.size());

        //then add per path variable
        for (size_t i=0;i<CurrSeparatrices.size();i++)
        {
            //get the two node indexes
            size_t Node0Index=CurrSeparatrices[i].first;
            size_t Node1Index=CurrSeparatrices[i].second;

            //safety checks
            assert(Node0Index>=0);
            assert(Node1Index>=0);
            assert(Node0Index<CoeffsV.size());
            assert(Node1Index<CoeffsV.size());

            //set the value of the constraint
            CoeffsV[Node0Index][i]=1;
            CoeffsV[Node1Index][i]=1;
        }
    }

    //set the conflict constraints
    void FindConflictsConstraints(std::vector<std::vector<ScalarType> > &CoeffsC,
                                  std::vector<ScalarType> &LhsC,
                                  std::vector<ScalarType> &RhsC)
    {
        CoeffsC.clear();
        LhsC.clear();
        RhsC.clear();

        //one row for each conflict, one column for each path variable
        CoeffsC.resize(CurrConflicts.size(),std::vector<ScalarType>(CurrPaths.size(),0));
        //set the boundary values
        LhsC.resize(CurrConflicts.size(),0);
        RhsC.resize(CurrConflicts.size(),1);

        for (size_t i=0;i<CurrConflicts.size();i++)
        {
            assert(CurrConflicts[i].Paths.size()>=2);
            for (size_t j=0;j<CurrConflicts[i].Paths.size();j++)
            {
                int IConfl=i;
                int IndexPath=CurrConflicts[i].Paths[j];
                assert(IndexPath>=0);
                assert(IndexPath<CurrPaths.size());
                //then set the constraint
                CoeffsC[IConfl][IndexPath]=1;
            }
        }
    }

    //update the solution and filter out the resolved ones
    void UpdateSolution(const std::vector<ScalarType> &solution,
                        std::vector<size_t> &UnsolvedNodes)
    {
        //       //clear the vector of non solved singularities
        //       CurrSingNodes.clear();

        UnsolvedNodes.clear();

        //the set of solved nodes
        std::set<int> SolvedVar;
        std::set<int> SolvedNodes;

        assert(solution.size()==CurrPaths.size());
        assert(solution.size()==CurrSeparatrices.size());

        //then update the solution
        for (size_t i=0;i<solution.size();i++)
        {
            int SolV=floor(solution[i]+0.5);
            if (SolV!=1)continue;

            //std::cout << "S " << solution[i] << std::endl;

            //add the current solution to the total solution
            TotalSolution.push_back(CurrPaths[i]);

            //then update the solved Nodes
            size_t NodeLocI0=CurrSeparatrices[i].first;
            size_t NodeLocI1=CurrSeparatrices[i].second;

            //then get the value of the Node in the Graph
            assert(NodeLocI0>=0);
            assert(NodeLocI0<CurrSingNodes.size());
            assert(NodeLocI1>=0);
            assert(NodeLocI1<CurrSingNodes.size());

            size_t NodeI0=CurrSingNodes[NodeLocI0];
            size_t NodeI1=CurrSingNodes[NodeLocI1];
            assert(NodeI0!=NodeI1);
            //should be found only once because of valence constraints
            assert(SolvedVar.count(NodeLocI0)==0);
            assert(SolvedVar.count(NodeLocI1)==0);
            assert(SolvedNodes.count(NodeI0)==0);
            assert(SolvedNodes.count(NodeI1)==0);

            //then insert the solved ones
            SolvedVar.insert(NodeLocI0);
            SolvedVar.insert(NodeLocI1);

            SolvedNodes.insert(NodeI0);
            SolvedNodes.insert(NodeI1);
        }

        //finally update the unsolved nodes
        for (size_t i=0;i<CurrSingNodes.size();i++)
        {
            size_t currINode=CurrSingNodes[i];
            if (SolvedNodes.count(currINode)>0)continue; //solved
            UnsolvedNodes.push_back(currINode);
        }
    }


    bool SolveILP(bool ConstrainSolution,
                  SepWeightType SepWType,
                  std::vector<ScalarType> &solution,
                  ScalarType &residual)
    {
        typedef SCIP_Solver<ScalarType> SolverType;
        typedef typename SolverType::VarType VarType;

        int t0 = clock();


        //create the objective function
        std::vector<ScalarType> ObjFunct;
        FindObjectiveFunct(ConstrainSolution,SepWType,ObjFunct);

        //the type of variables
        std::vector<VarType> varTypes(CurrPaths.size(),VarType::BINARY);

        std::vector<ScalarType> lowerBounds(CurrPaths.size(),0);
        std::vector<ScalarType> upperBounds(CurrPaths.size(),1);

        //constraint for singularity valence
        std::vector<std::vector<ScalarType> > CoeffsV;
        std::vector<ScalarType> LhsV,RhsV;
        FindSingularityConstraints(CoeffsV,LhsV,RhsV,ConstrainSolution);

        //constraints for conflicts
        std::vector<std::vector<ScalarType> > CoeffsC;
        std::vector<ScalarType> LhsC;
        std::vector<ScalarType> RhsC;
        std::cout << "Finding Conflicts " << std::endl;
        fflush(stdout);
        FindConflictsConstraints(CoeffsC,LhsC,RhsC);

        //then instantiate the problem
        typedef SCIP_Solver<ScalarType> SolverType;
        typedef typename SolverType::VarType VarType;

        //create the solver
        SolverType Solver;

        std::cout << "Creating Problem" << std::endl;
        fflush(stdout);

        if (ConstrainSolution)
            Solver.newProblem(ObjFunct,varTypes,lowerBounds,upperBounds,false);
        else
            Solver.newProblem(ObjFunct,varTypes,lowerBounds,upperBounds,true);

        //add the constraints
        std::cout << "Adding Constraints" << std::endl;
        fflush(stdout);
        std::cout << "Add 1" << std::endl;

        Solver.addConstraints(CoeffsV,LhsV,RhsV);

        std::cout << "Add 2" << std::endl;
        Solver.addConstraints(CoeffsC,LhsC,RhsC);

        int t1 = clock();

        std::cout << "Global Solve Step " << std::endl;
        std::cout << "Num Variables " << CurrPaths.size() << std::endl;
        std::cout << "Num Constraints " << CoeffsV.size()+CoeffsC.size() << std::endl;
        std::cout << "Setup Time " << (float)(t1 - t0)/CLOCKS_PER_SEC << "sec" << std::endl;
        fflush(stdout);

        bool solved=Solver.solve(solution,residual);
        int t2 = clock();

        if (solved) std::cout << "Solved Successfully !! " << std::endl;
        else std::cout << "Not solved !! " << std::endl;

        std::cout << "Solve Time " << (float)(t2 - t1)/CLOCKS_PER_SEC << "sec" << std::endl;

        return solved;
    }


    bool SolveStep(size_t maxConn,
                   ScalarType MaxDev,
                   ScalarType PenaltyDrift,
                   SepWeightType SepWType,
                   bool progressivePaths,
                   bool ConstrainSolution)
    {
        //initialize the new separatrix
        bool found=FindSeparatrixStep(maxConn,MaxDev,PenaltyDrift,progressivePaths);
        if (!found)return false;

        //then update the conflicts
        ConflictFinder<MeshType> ConflFinder(anigraph);
        ConflFinder.FindConflicts(CurrPaths,CurrConflicts,true);

        //then perform the solving step
        std::vector<ScalarType> solution;
        ScalarType residual;

        if (ConstrainSolution)std::cout << "Solving Constrained " << std::endl;
        if (!ConstrainSolution)std::cout << "Solving Unconstrained " << std::endl;
        bool solved=SolveILP(ConstrainSolution,SepWType,solution,residual);

        if (!solved)return false;//not found a feasible solution

        size_t NumUnSolved0=CurrSingNodes.size();

        //and update the solution
        std::vector<size_t> UnsolvedNodes;
        UpdateSolution(solution,UnsolvedNodes);

        std::cout << "Unsolved " << UnsolvedNodes.size() << " Singularity Nodes" << std::endl;

        if (!progressivePaths)
        {
            CurrSingNodes=std::vector<size_t>(UnsolvedNodes.begin(),UnsolvedNodes.end());
            size_t NumUnSolved1=CurrSingNodes.size();
            assert(NumUnSolved1<=NumUnSolved0);
            //write some stats

            //then set the barrier for the next step
            SetSolutionAsBarrier(false);
            return (NumUnSolved1<NumUnSolved0);
        }else
        {
            if (UnsolvedNodes.size()==0)return false;
            const size_t MaxMult=4;
            bool incremented=false;
            for (size_t j=0;j<UnsolvedNodes.size();j++)
            {
                size_t currNode=UnsolvedNodes[j];
                assert(SingToIndex.count(currNode)>0);
                int curr_index=SingToIndex[currNode];
                size_t newVal=TargetEvaluated[curr_index]*2;
                if (newVal<=(maxConn*MaxMult))
                {
                    TargetEvaluated[curr_index]=newVal;
                    incremented=true;
                }
            }
        }
    }

    void UpdatateUnsolvedNodes()
    {
        UnsolvedNodes.clear();

        std::set<size_t> VisitedN;
        //the vector of visited sinks
        for (size_t i=0;i<TotalSolution.size();i++)
        {
            //get the two nodes
            size_t Node0=TotalSolution[i].nodes[0];
            size_t Node1=TotalSolution[i].nodes.back();

            //if it is not a singularity or a sink
            //then it is part of a t-junction
            if (!anigraph.IsEdgeNode(Node0))
            {
                if (anigraph.IsSink(Node0)) Node0=anigraph.OppositeSink(Node0);
                assert(anigraph.IsSingularity(Node0));
                VisitedN.insert(Node0);
            }
            if (!anigraph.IsEdgeNode(Node1))
            {
                if (anigraph.IsSink(Node1)) Node1=anigraph.OppositeSink(Node1);
                assert(anigraph.IsSingularity(Node1));
                VisitedN.insert(Node1);
            }
        }
        for (size_t i=0;i<InitSingNodes.size();i++)
        {
            size_t IndexSing=InitSingNodes[i];
            if (VisitedN.count(IndexSing)>0)continue;
            UnsolvedNodes.push_back(IndexSing);
        }
    }

    void PrintStats()
    {
        std::cout << "*** SOME STATISTIC ***" << std::endl;


        UpdatateUnsolvedNodes();
        std::cout << "There are " << UnsolvedNodes.size() << " Unsolved Nodes" << std::endl;

        //then print the conflicts
        ConflictFinder<MeshType> ConflFinder(anigraph);
        std::vector<ConflictInfo> TempConflicts;
        ConflFinder.FindConflicts(TotalSolution,TempConflicts,true);
        /*std::cout << "There are " << TempConflicts.size() << " Conflicts" << std::endl;*/

        std::cout << "*** END STATISTIC ***" << std::endl;

    }

    void AddTJunctionsStep(ScalarType MaxDev,
                           ScalarType PenaltyDrift,
                           int Ndir=2,
                           bool onlyBorder=false)
    {
        //first erase all blocks
        anigraph.SetAllValid();

        //block close to old solution
        SetSolutionAsBarrier(true);

        //deselect all nodes for visualization purposes
        anigraph.DeselectAllNodes();

        //collect all sources
        std::vector<size_t> SourcesN,BorderN;
        if (!onlyBorder)
        anigraph.TangentialNodes(TotalSolution,SourcesN,Ndir);

        //add border Nodes
        anigraph.GetBorderNodes(BorderN);
        SourcesN.insert(SourcesN.end(),BorderN.begin(),BorderN.end());

        //then call dyjkstra
        std::vector<PathType> TPaths;
        //anigraph.FindPatternToSingularities(SourcesN,MaxDev,PenaltyDrift,TPaths);
        AnisotropicQuery<MeshType>::FindPatternToSingularities(anigraph,SourcesN,MaxDev,PenaltyDrift,TPaths);

        //then invert each path
        for (size_t i=0;i<TPaths.size();i++)
        {
            assert(anigraph.IsSink(TPaths[i].nodes.back()));
            //invert the path
            anigraph.InvertPath(TPaths[i]);
        }

        //check for conflicts
        ConflictFinder<MeshType> ConflFinder(anigraph);
        std::vector<ConflictInfo> TempConflicts;
        ConflFinder.FindConflicts(TPaths,TempConflicts,true);
        //assert(TPaths.size()==TempConflicts.size());

        //create a set of all path that has at least one conflict
        struct ConflicPriority
        {
            size_t IndexPath;
            std::vector<size_t> Conflicting;

//            ConflicPriority(std::vector<size_t> &_Conflicting,size_t _IndexPath)
//            {
//                IndexPath=_IndexPath;
//                Conflicting=std::vector<size_t>(_Conflicting.begin(),_Conflicting.end());
//            }


            ConflicPriority(size_t _IndexPath)
            {
                IndexPath=_IndexPath;
            }

            inline bool operator <(const ConflicPriority &Cp)const
            {return (Conflicting.size()<Cp.Conflicting.size());}
        };

        std::vector<ConflicPriority> ConflictStack;
        for (size_t i=0;i<TPaths.size();i++)
            ConflictStack.push_back(ConflicPriority(i));

        std::vector<std::pair<size_t,size_t> > ConflPairs;
        for (int i=0;i<TempConflicts.size();i++)
        {
            //for each conflicting pair insert into the structure
           for (size_t j=0;j<TempConflicts[i].Paths.size()-1;j++)
               for (size_t k=j+1;k<TempConflicts[i].Paths.size();k++)
               {
                  size_t ConflP0=TempConflicts[i].Paths[j];
                  size_t ConflP1=TempConflicts[i].Paths[k];
                  assert(ConflP0!=ConflP1);
                  ConflPairs.push_back(std::pair<size_t,size_t> (std::min(ConflP0,ConflP1),std::max(ConflP0,ConflP1)));
               }
        }
        std::sort(ConflPairs.begin(),ConflPairs.end());
        std::vector<std::pair<size_t,size_t> >::iterator last= std::unique(ConflPairs.begin(), ConflPairs.end());
        ConflPairs.erase(last, ConflPairs.end());

        //then insert into the vector
        for (int i=0;i<ConflPairs.size();i++)
        {
            ConflictStack[ConflPairs[i].first].Conflicting.push_back(ConflPairs[i].second);
            ConflictStack[ConflPairs[i].second].Conflicting.push_back(ConflPairs[i].first);
        }

        assert(ConflictStack.size()==TPaths.size());

        std::sort(ConflictStack.begin(),ConflictStack.end());


        std::vector<int> ToInsert(TPaths.size(),-1);
        for (size_t i=0;i<ConflictStack.size();i++)
        {
            int IndexP0=ConflictStack[i].IndexPath;
            if (ToInsert[IndexP0]==0)continue;//already disabled for some other path
            assert(ToInsert[IndexP0]==-1);

            //then insert it
            ToInsert[IndexP0]=1;

            //and disable inserting from the rest
            for (size_t j=0;j<ConflictStack[i].Conflicting.size();j++)
            {
                int IndexP1=ConflictStack[i].Conflicting[j];
                ToInsert[IndexP1]=0;
            }
        }

        //then add one by one the one that are suitable
        for (int i=0;i<TPaths.size();i++)
        {
            assert(ToInsert[i]!=-1);
            if (ToInsert[i]==0)continue;
            TotalSolution.push_back(TPaths[i]);
            TotalSolution.back().isTJunction=true;
        }
//        //no conflicts or at least 2 paths
//        assert((ConflSet.size()==0)||(ConflSet.size()>=2));

//        //then add the ones that have no conflict and only one that have conflicts
//        // to go on for the next cycle
//        for (int i=0;i<TPaths.size();i++)
//        {
//            if (ConflSet.count(i)>0)continue;
//            //and add to the solution
//            TotalSolution.push_back(TPaths[i]);
//            TotalSolution.back().isTJunction=true;
//        }

//        if (TempConflicts.size()==0)return;

//        size_t firstIndex=TempConflicts[0].Paths[0];

//        assert(firstIndex>=0);
//        assert(firstIndex<TPaths.size());

//        TotalSolution.push_back(TPaths[firstIndex]);

        //return true;
    }

//    void AddTJunctionsFrom(int IndexNode,
//                           ScalarType MaxDev,
//                           ScalarType PenaltyDrift,
//                           int Ndir=2,
//                           bool onlyBorder=false)
//    {
//        //first erase all blocks
//        anigraph.SetAllValid();

//        //block close to old solution
//        SetSolutionAsBarrier(true);

//        //deselect all nodes for visualization purposes
//        anigraph.DeselectAllNodes();

//        //collect all targets
//        std::vector<size_t> TargetN,BorderN;
//        if (!onlyBorder)
//        anigraph.TangentialNodes(TotalSolution,TargetN,Ndir);

//        //add border Nodes
//        anigraph.GetBorderNodes(BorderN);
//        SourcesN.insert(SourcesN.end(),BorderN.begin(),BorderN.end());

//        //then call dyjkstra
//        std::vector<PathType> TPaths;
//        //anigraph.FindPatternToSingularities(SourcesN,MaxDev,PenaltyDrift,TPaths);
//        AnisotropicQuery<MeshType>::FindPatternToSingularities(anigraph,SourcesN,MaxDev,PenaltyDrift,TPaths);
//        AnisotropicQuery<MeshType>::computePerNodeDijsktra(anigraph,SourcesN,MaxDev,PenaltyDrift,TPaths);

//        //then invert each path
//        for (size_t i=0;i<TPaths.size();i++)
//        {
//            assert(anigraph.IsSink(TPaths[i].nodes.back()));
//            //invert the path
//            anigraph.InvertPath(TPaths[i]);
//        }

//        //check for conflicts
//        ConflictFinder<MeshType> ConflFinder(anigraph);
//        std::vector<ConflictInfo> TempConflicts;
//        ConflFinder.FindConflicts(TPaths,TempConflicts,true);

//        //create a set of all path that has at least one conflict
//        std::set<size_t> ConflSet;
//        for (int i=0;i<TempConflicts.size();i++)
//            for (int j=0;j<TempConflicts[i].Paths.size();j++)
//                ConflSet.insert(TempConflicts[i].Paths[j]);

//        //no conflicts or at least 2 paths
//        assert((ConflSet.size()==0)||(ConflSet.size()>=2));

//        //then add the ones that have no conflict and only one that have conflicts
//        // to go on for the next cycle
//        for (int i=0;i<TPaths.size();i++)
//        {
//            if (ConflSet.count(i)>0)continue;
//            //and add to the solution
//            TotalSolution.push_back(TPaths[i]);
//            TotalSolution.back().isTJunction=true;
//        }
//        if (TempConflicts.size()==0)return;

//        size_t firstIndex=TempConflicts[0].Paths[0];

//        assert(firstIndex>=0);
//        assert(firstIndex<TPaths.size());

//        TotalSolution.push_back(TPaths[firstIndex]);

//        //return true;
//    }

    void EsteemCutOff(int NumNeigh,
                      ScalarType MaxDev,
                      ScalarType PenaltyDrift)
    {
        anigraph.SetAllValid();

        std::cout << "*** ESTEEEM CUTOFF DISTANCE ***" << std::endl;
        const int MaxTest=10;
        //get the singularities and sinks informations
        std::vector<std::vector<size_t > > SingNodes,SinkNodes;
        std::vector<size_t> SingVertex,ExpValence;
        anigraph.GetSingularitiesInfo(SingNodes,SinkNodes,SingVertex,ExpValence);

        std::vector<size_t> testSing;
        //then initialize the vector of initial singularities
        for (size_t i=0;i<SingNodes.size();i++)
            for (size_t j=0;j<SingNodes[i].size();j++)
                testSing.push_back(SingNodes[i][j]);

        //then for each singularity find the possible connections to sinks
        cutoffL=0;
        int Num=0;
        for (size_t i=0;i<testSing.size();i++)
        {
            //associate the index to the node
            size_t NodeSing=testSing[i];
            //SingToIndex[NodeSing]=i;

            //safety check that is a singularity
            assert(anigraph.IsSingularity(NodeSing));
            std::vector<PathType> TempPaths;
            AnisotropicQuery<MeshType>::FindPatternToSingularities(anigraph,NodeSing,MaxDev,PenaltyDrift,TempPaths,MaxTest,false,std::numeric_limits<ScalarType>::max());
            //anigraph.FindPatternToSingularities(NodeSing,MaxDev,PenaltyDrift,TempPaths,MaxTest,false,std::numeric_limits<ScalarType>::max());
            for (size_t j=0;j<TempPaths.size();j++)
                cutoffL+=TempPaths[j].length;
            Num+=TempPaths.size();
        }
        cutoffL/=(ScalarType)Num;
        cutoffL*=((ScalarType)NumNeigh)/((ScalarType)MaxTest);
        std::cout << "*** ESTEEEMATED" << cutoffL <<" ***" << std::endl;
    }



    void UpdateCandidatesEdgeM()
    {
        CandidateEMesh.Clear();

        //get the singularities
        std::vector<size_t > SingNode,SinkNode;
        anigraph.GetSingularityNodes(SingNode,SinkNode);

        //then create the edgemesh
        for (size_t i=0;i<AllCandidates.size();i++)
        //for (size_t i=0;i<TotalSolution.size();i++)
        {
            bool to_add=true;
            if (OldSingCandidate!=-1)//in this case check the singularity index
            {
                int Index=(OldSingCandidate % SingNode.size());
                int TargetSing=SingNode[Index];
                size_t SingNode0,SingNode1;
                anigraph.GetSingExtremes(AllCandidates[i],SingNode0,SingNode1);
                if ((SingNode0!=TargetSing)&&(SingNode1!=TargetSing))to_add=false;
            }
            if (!to_add)continue;
            std::vector<PathType> currP(1,AllCandidates[i]);

            MyEMesh swapM;
            std::vector<bool> IsLoop(currP.size(),false);
            anigraph.GetEdgeMesh(currP,IsLoop,swapM,anigraph.Mesh());
            for (size_t j=0;j<swapM.edge.size();j++)
                swapM.edge[j].PathIndex=i;
            vcg::tri::Append<MyEMesh,MyEMesh>::Mesh(CandidateEMesh,swapM);
        }

        //
        CandidateEMesh.SelectSingularVert();
        //then smooth locally
        CandidateEMesh.InitFaceEdgeMap(anigraph.Mesh());
        CandidateEMesh.ComputeElevation(anigraph.Mesh());

        CandidateEMesh.SmoothLocally(anigraph.Mesh(),20,false);
        CandidateEMesh.SelectSingularVert();

        CandidateEMesh.InitColorByPath();

        //then elevate
        CandidateEMesh.ApplyElevation();
        CandidateEMesh.Smooth(5);
    }

    void  AddTJunctionsCycle(ScalarType MaxDev,
                        ScalarType PenaltyDrift,
                        int Ndir=2,
                        bool only_border=false)
    {
        int t0=clock();
        UpdatateUnsolvedNodes();
        size_t unsolved0=UnsolvedNodes.size();
        size_t unsolved1=unsolved0;
        do{
            unsolved0=unsolved1;
            AddTJunctionsStep(MaxDev,PenaltyDrift,2,only_border);
            UpdatateUnsolvedNodes();
            unsolved1=UnsolvedNodes.size();
            std::cout << "*** TSTEP ***" << std::endl;
        }while (unsolved1<unsolved0);
        PrintStats();
        Stats.UnresolvedDirections1=UnsolvedNodes.size();
        int t1=clock();
        Stats.timeSeparatrix+=(ScalarType)(t1-t0)/(ScalarType)(CLOCKS_PER_SEC);
    }

public:

    int MaxIterations;
    int OldSingCandidate;

    void  AddTJunctions(ScalarType MaxDev,
                        ScalarType PenaltyDrift,
                        int Ndir=2)
    {
        AddTJunctionsCycle(MaxDev,PenaltyDrift,Ndir,true);
        AddTJunctionsCycle(MaxDev,PenaltyDrift,Ndir,false);
    }

    void Reset()
    {
        //the first time must be initialized
        OldSingCandidate=-2;
        InitSingNodes.clear();
        TotalSolution.clear();
        CurrSingNodes.clear();
        CurrSeparatrices.clear();
        CurrPaths.clear();
        CurrConflicts.clear();
        anigraph.SetAllValid();
        NumEvaluated.clear();
        TargetEvaluated.clear();
        AllCandidates.clear();

        minL=std::numeric_limits<ScalarType>::max();
        maxL=0;
    }

#ifndef NO_TRACING_OPENGL

    void GLDrawCandidates(int curr_candidate)
    {
        if (curr_candidate!=OldSingCandidate)
        {
            OldSingCandidate=curr_candidate;
            UpdateCandidatesEdgeM();
        }
        CandidateEMesh.GLDraw();
    }

    //draw a given path,
    //the singularities it is involved and its conflicts
    void DrawSolution(bool show_only_t_junctions=false)
    {
        if (TotalSolution.size()==0)return;
        size_t sizeSol=TotalSolution.size();
        for (int i=0;i<sizeSol;i++)
        {
            vcg::Color4b currC=vcg::Color4b::Scatter(sizeSol,i);

            if ((show_only_t_junctions)&&(!TotalSolution[i].isTJunction))continue;

            assert(TotalSolution[i].nodes.size()>1);
            anigraph.GlDrawPath(TotalSolution[i],currC,20);
        }
    }


    void GLDrawUnsolved(ScalarType MaxDev,
                        ScalarType PenaltyDrift)
    {
        for (size_t i=0;i<UnsolvedNodes.size();i++)
        {
            int IndexN=UnsolvedNodes[i];
            anigraph.GlDrawNodeDebugInfo(IndexN,MaxDev,PenaltyDrift,0.001,false,true);
        }
    }
#endif

    void SolveGlobal(size_t maxConn,
                     ScalarType MaxDev,
                     ScalarType PenaltyDrift,
                     SepWeightType SepWType,
                     SolveMode SolvMode)
    {

        cutoffL=-1;
        if (PercentileCutOff)EsteemCutOff(maxConn,MaxDev,PenaltyDrift);

        Stats.NumNeighbors=maxConn;

        ScalarType t0=clock();
        Reset();

        Initialize();

        int iteration=0;
        //check if it is constrained than do a single step
        switch ( SolvMode )
        {
        case SMConstraintSol:
        {
            iteration++;
            SolveStep(maxConn,MaxDev,PenaltyDrift,SepWType,false,true);
            break;
        }
        case SMIterativeConstr:
        {
            bool progressed=true;
            do
            {
                iteration++;
                progressed&=SolveStep(maxConn,MaxDev,PenaltyDrift,SepWType,false,false);
                maxConn*=4;
                cutoffL*=4;
            }while ((progressed)&&(iteration<MaxIterations));
            break;
        }
        case SMItaretiveGlobal:
        {
            //then initialize the target and the inital sampled
            NumEvaluated.resize(CurrSingNodes.size(),0);
            TargetEvaluated.resize(CurrSingNodes.size(),maxConn);
            bool progressed=true;
            do
            {
                iteration++;
                progressed&=SolveStep(maxConn,MaxDev,PenaltyDrift,SepWType,true,false);
                if (!progressed)
                    SolveStep(maxConn,MaxDev,PenaltyDrift,SepWType,false,true);
            }while ((progressed)&&(iteration<MaxIterations));
            break;
        }
        }
        PrintStats();

        Stats.UnresolvedDirections0=UnsolvedNodes.size();
        Stats.NumIterations=iteration;
        ScalarType t1=clock();
        Stats.timeSeparatrix=(ScalarType)(t1-t0)/(ScalarType)(CLOCKS_PER_SEC);

        //UpdateCandidatesEdgeM();
    }

};
#endif
