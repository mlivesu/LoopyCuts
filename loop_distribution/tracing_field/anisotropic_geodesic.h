#ifndef ANISOTROPIC_GEODESIC_QUERY
#define ANISOTROPIC_GEODESIC_QUERY

#include "graph/anisotropic_graph.h"
#include "loop_common_functions.h"

template < class MeshType >
class AnisotropicQuery
{
    typedef AnisotropicGraph<MeshType> AnisotropicGraph;

    typedef typename MeshType::CoordType    CoordType;
    typedef typename MeshType::VertexType   VertexType;
    typedef typename MeshType::FaceType     FaceType;
    typedef typename MeshType::ScalarType   ScalarType;

    //typedef typename BaseGeodesicNode::NeighborsIterator NeighborsIterator;

    /*typedef EdgeGeodesicNode<CoordType> EdgeGeodesicNode;
    typedef SingularityGeodesicNode<CoordType> SingularityGeodesicNode;
    typedef SinkGeodesicNode<CoordType> SinkGeodesicNode;

    typedef NeighborInfo<ScalarType> NeighborInfo;
    typedef BaseGeodesicNode<CoordType> BaseGeodesicNode;*/

    //public:
    typedef Path<ScalarType> Path;
    typedef BaseGeodesicNode<CoordType> BaseGeodesicNode;
    typedef typename BaseGeodesicNode::NeighborsIterator NeighborsIterator;



    /**
     * @brief The NodeDist struct
     */
    struct NodeDist {
        size_t                      NIndex; ///< Node index.
        std::vector<ScalarType>    *Dist;   ///< Distance of every node from a given seed.

        /**
         * @brief NodeDist Constructor.
         * @param _NIndex
         * @param _Dist
         */
        NodeDist(const size_t _NIndex,
                 std::vector<ScalarType> *_Dist)
            : NIndex(_NIndex),
              Dist(_Dist) { }

        /**
         * @brief NodeDist Copy constructor.
         * @param nd
         */
        NodeDist(const NodeDist &nd)
            : NIndex(nd.NIndex),
              Dist(nd.Dist) { }

        /**
         * @brief operator = Copy assignment operator.
         * @param nd
         * @return
         */
        NodeDist & operator =(const NodeDist &nd) {
            if (this != &nd) {
                NIndex = nd.NIndex;
                Dist = nd.Dist;
            }
            return *this;
        }

        /**
         * @brief operator <
         * @param o
         * @return
         */
        bool operator < (const NodeDist &o) const {
            if( (*Dist)[NIndex] != (*o.Dist)[o.NIndex])
                return (*Dist)[NIndex] > (*o.Dist)[o.NIndex];
            return NIndex < o.NIndex;
        }
    };

public:

    struct DijkstraParam
    {
        std::vector<std::vector<size_t> > sourceGroups;
        ScalarType maxDist;
        size_t target;
        size_t maxSink;

        DijkstraParam()
        {
            maxDist=std::numeric_limits<ScalarType>::max();
            target=std::numeric_limits<size_t>::max();
            maxSink=std::numeric_limits<size_t>::max();
        }
    };

    static void retrievePath(AnisotropicGraph &AniGraph,
                             const size_t startNode,
                             const size_t endNode,
                             std::vector<size_t> &path)
    {
        // get the nodes from
        size_t curr_node = endNode;
        size_t prev_node = curr_node;
        std::list<size_t> nodeList;
        do {
            nodeList.push_back(curr_node);
            prev_node = curr_node;
            curr_node = AniGraph.Father[curr_node];
        } while (curr_node != std::numeric_limits<size_t>::max() &&
                 curr_node != startNode &&
                 curr_node != endNode &&
                 curr_node != prev_node);
        if ((curr_node == startNode)&&
                (startNode!=std::numeric_limits<size_t>::max()))
            nodeList.push_back(startNode);
        path.assign(nodeList.rbegin(), nodeList.rend());
    }

    static void retrievePath(AnisotropicGraph &AniGraph,
                             const size_t endNode,
                             std::vector<size_t> &path)
    {
        size_t startNode=std::numeric_limits<size_t>::max();
        retrievePath(AniGraph,startNode,endNode,path);
    }

    static void retrievePath(AnisotropicGraph &AniGraph,
                             const size_t endNode,
                             Path &retPath)
    {
        retPath.length = AniGraph.Distance[endNode];
        retrievePath(AniGraph,endNode,retPath.nodes);
    }

    /**
     * @brief computePerNodeDijsktra
     * @param seedNodes
     * @param maxDistanceThr
     * @param target
     * @param max_deviation
     * @param penalty_drift
     */
    static int computePerNodeDijsktra(AnisotropicGraph &AniGraph,
                                      const DijkstraParam &Param,
                                      const ScalarType max_deviation = ScalarType(45),
                                      const ScalarType penalty_drift = ScalarType(2),
                                      std::vector<size_t> *withinRadius=NULL,
                                      std::vector<ScalarType> *PrevDistance=NULL,
                                      std::vector<size_t> *BoderNode=NULL)
    {
        // initialize the vector of fathers
        //size_t t0=clock();
        AniGraph.Group.assign(AniGraph.Graph.size(), std::numeric_limits<size_t>::max());
        AniGraph.Father.assign(AniGraph.Graph.size(), std::numeric_limits<size_t>::max());
        AniGraph.Source.assign(AniGraph.Graph.size(), std::numeric_limits<size_t>::max());
        //size_t t1=clock();
        //std::cout<<"Init Structure Time "<<(t1-t0)/(ScalarType)CLOCKS_PER_SEC<<std::endl;

        if (PrevDistance!=NULL)
        {
            assert(PrevDistance->size()==AniGraph.Graph.size());
            AniGraph.Distance.assign(PrevDistance->begin(),PrevDistance->end());
        }
        else
            AniGraph.Distance.assign(AniGraph.Graph.size(), std::numeric_limits<ScalarType>::max());

        std::deque<NodeDist> heap;

        // push all seeds in the heap
        for (size_t i = 0; i < Param.sourceGroups.size(); ++i)
            for (size_t j = 0; j < Param.sourceGroups[i].size(); ++j)
            {
                size_t NodeI=Param.sourceGroups[i][j];
                assert(AniGraph.Graph[NodeI]->Valid);
                AniGraph.Distance[NodeI] = ScalarType(0);
                AniGraph.Father[NodeI] = NodeI;
                AniGraph.Source[NodeI] = NodeI;
                AniGraph.Group[NodeI]=i;
                heap.push_back(NodeDist(NodeI, &AniGraph.Distance));
                if (withinRadius==NULL)continue;
                withinRadius->push_back(NodeI);
            }

        std::make_heap(heap.begin(), heap.end());

        std::set<size_t> FoundSink;

        // explore the graph
        while (!heap.empty())
        {
            // pop from the heap
            std::pop_heap(heap.begin(), heap.end());
            size_t currNode = heap.back().NIndex;

            //if not initialized
            if (AniGraph.IsSink(currNode)&&
                    (Param.maxSink!=std::numeric_limits<size_t>::max()))
                FoundSink.insert(currNode);

            if (FoundSink.size()==Param.maxSink)return currNode;

            heap.pop_back();

            // finish when reach the target
            if (Param.target != std::numeric_limits<size_t>::max() &&
                    Param.target == currNode &&
                    AniGraph.Distance[currNode] > ScalarType(0))
                return currNode;

            //check gruop and opposites
            if (AniGraph.IsEdgeNode(currNode)&&(!AniGraph.IsBorderNode(currNode)))
            {
                size_t CurrGroup=AniGraph.Group[currNode];
                int OppNode=AniGraph.OppositeNode(currNode);
                assert(OppNode!=-1);
                size_t OppGroup=AniGraph.Group[OppNode];

                //two different fronts meets!
                if ((CurrGroup!=std::numeric_limits<size_t>::max())&&
                        (OppGroup!=std::numeric_limits<size_t>::max())&&
                        (AniGraph.Distance[currNode] > ScalarType(0))&&
                        (AniGraph.Distance[OppNode] > ScalarType(0))&&
                        (CurrGroup!=OppGroup))
                    return currNode;
            }



            // look at the neighbors
            for(NeighborsIterator it = AniGraph.Graph[currNode]->neighbors.begin(); it != AniGraph.Graph[currNode]->neighbors.end(); ++it)
            {

                //discard if it is not valid
                if (!it->Valid)
                {
                    if (BoderNode!=NULL)
                        BoderNode->push_back(currNode);
                    continue;
                }


                int IndexNeigh=it->nodeID;
                assert(IndexNeigh < (int)AniGraph.Graph.size());


                if (it->angle > max_deviation)
                    continue;


                // compute the new distance to the current neighbor
                ScalarType nextDistance = AniGraph.Distance[currNode] + it->geoDistance(max_deviation, penalty_drift);

                // discard neighbors too much far
                if (nextDistance >= Param.maxDist)
                    continue;

                if (withinRadius!=NULL)
                    withinRadius->push_back(currNode);

                assert(AniGraph.Graph[IndexNeigh]->Valid);

                // update the new distance and push in the heap
                //if (Distance[IndexNeigh] <= ScalarType(0) || nextDistance < Distance[IndexNeigh])
                if (nextDistance < AniGraph.Distance[IndexNeigh])
                {
                    AniGraph.Father[IndexNeigh] = currNode;
                    AniGraph.Source[IndexNeigh] = AniGraph.Source[currNode];
                    AniGraph.Distance[IndexNeigh] = nextDistance;
                    AniGraph.Group[IndexNeigh]=AniGraph.Group[AniGraph.Source[IndexNeigh]];
                    assert((int)AniGraph.Group[IndexNeigh]!=-1);

                    heap.push_back(NodeDist(IndexNeigh, &AniGraph.Distance));
                    std::push_heap(heap.begin(), heap.end());
                }
            }
        }
        return -1;

    }



    static void FindNodesCloseToPath(AnisotropicGraph &AniGraph,
                                     const Path &path,
                                     const ScalarType radius,
                                     const ScalarType max_deviation,
                                     const ScalarType penalty_drift,
                                     std::vector<size_t> &closeNodes,
                                     bool bidirectional=false)
    {
        closeNodes.clear();
        DijkstraParam DParam;
        DParam.maxDist=radius;
        DParam.sourceGroups.push_back(std::vector<size_t>());
        for (size_t i=0;i<path.nodes.size();i++)
        {
            DParam.sourceGroups[0].push_back(path.nodes[i]);
            //            if (!bidirectional)continue;
            //            int OppI=OppositeNode(path.nodes[i]);
            //            DParam.sourceGroups[0].push_back(OppI);
        }
        std::vector<size_t> closeNodesDir0;
        computePerNodeDijsktra(AniGraph,DParam,max_deviation,penalty_drift,&closeNodesDir0);
        //then add the opposites
        for (size_t i=0;i<closeNodesDir0.size();i++)
        {
            int IndexN=closeNodesDir0[i];
            closeNodes.push_back(IndexN);
            if (!bidirectional)continue;
            if (AniGraph.IsSingularity(IndexN)||AniGraph.IsSink(IndexN))continue;

            //add the opposite on the same face
            int OppSameF=AniGraph.OppositeNodeSameF(IndexN);
            if (OppSameF!=-1)closeNodes.push_back(OppSameF);

            //get the opposite
            int OppI=AniGraph.OppositeNode(closeNodesDir0[i]);
            closeNodes.push_back(OppI);
            OppSameF=AniGraph.OppositeNodeSameF(OppI);
            if (OppSameF!=-1)closeNodes.push_back(OppSameF);
        }
        std::sort(closeNodes.begin(),closeNodes.end());
        closeNodes.erase( std::unique( closeNodes.begin(), closeNodes.end() ), closeNodes.end() );

    }

    /**
     * @brief FindPatternToSingularities finds at most maxNPaths paths from startNode to each singularity.
     * @param startNode
     * @param max_deviation
     * @param penalty_drift
     * @param paths
     * @param maxNPaths
     */
    static void FindPatternToSingularities(AnisotropicGraph &AniGraph,
                                           const size_t startNode,
                                           const ScalarType max_deviation,
                                           const ScalarType penalty_drift,
                                           std::vector<Path> &paths,
                                           const size_t maxNPaths,
                                           bool printmsg,
                                           ScalarType MaxDist)
    {
        // compute distances
        if (printmsg)
        {
            std::cout << "Computing Dijkstra ... ";
            std::cout.flush();
        }
        int t0 = clock();
        DijkstraParam Param;
        Param.sourceGroups.resize(1);
        Param.sourceGroups[0].push_back(startNode);
        Param.maxSink=maxNPaths;
        Param.maxDist=MaxDist;
        //std::vector<std::vector<size_t> seeds(1, startNode);
        computePerNodeDijsktra(AniGraph,Param,max_deviation, penalty_drift);
        int t1 = clock();

        if (printmsg)
        {
            std::cout << "done. Time " << (float)(t1 - t0)/CLOCKS_PER_SEC << std::endl;

            // extract reverse path
            std::cout << ". Extracting paths to sinks ... ";
            std::cout.flush();
        }
        t0 = clock();
        //        std::map<CoordType, Path> perSinkPaths;
        std::deque<Path> nearestPaths;
        typename std::deque<Path>::iterator lowerIt;
        Path newPath;
        for (size_t i = AniGraph.nEdgeNodes + AniGraph.nSingNodes; i < AniGraph.Graph.size(); ++i)
        {
            assert(AniGraph.Graph[i]->nodeType == SinkNodeType);
            if (!AniGraph.Graph[i]->Valid)continue;
            newPath.clear();
            newPath.length = AniGraph.Distance[i];
            retrievePath(AniGraph,startNode, i, newPath.nodes);
            if(newPath.nodes.size()<=1){
                //                printf("WARNING, NON REACHEABLE SINK\n");
                //                fflush(stdout);
                continue;
            }
            //check if it start and end in the same node

            lowerIt = std::lower_bound(nearestPaths.begin(), nearestPaths.end(), newPath);
#if __cplusplus >= 201103L
            nearestPaths.insert(lowerIt, std::move(newPath));
#else
            nearestPaths.insert(lowerIt, newPath);
#endif
        }

        // copy out the paths
        paths.resize(std::min(nearestPaths.size(), maxNPaths));
        for (size_t i = 0; i < paths.size(); ++i)
#if __cplusplus >= 201103L
            paths[i] = std::move(nearestPaths[i]);
#else
            paths[i] = nearestPaths[i];
#endif
        t1 = clock();
        if (printmsg)
            std::cout << "done. Time " << (float)(t1 - t0)/CLOCKS_PER_SEC << std::endl;
    }

    static bool FindShortestPatternBetween(AnisotropicGraph &AniGraph,
                                           std::vector<size_t> startNodes,
                                           std::vector<size_t> TargetNodes,
                                           const ScalarType max_deviation,
                                           const ScalarType penalty_drift,
                                           Path &CurrPath)
    {
        DijkstraParam Param;
        Param.sourceGroups.resize(2);
        Param.sourceGroups[0]=startNodes;
        Param.sourceGroups[1]=TargetNodes;
        int middleN=computePerNodeDijsktra(AniGraph,Param,max_deviation, penalty_drift);
        if (middleN==-1)return false;

        CurrPath.nodes.clear();

        //get the opposite node
        int OppMiddle=AniGraph.OppositeNode(middleN);
        assert(OppMiddle!=-1);

        //set the lenght
        CurrPath.length=AniGraph.Distance[middleN];
        CurrPath.length+=AniGraph.Distance[OppMiddle];

        //check group consistency
        assert(AniGraph.Group[middleN]!=AniGraph.Group[OppMiddle]);
        assert((AniGraph.Group[middleN]==0)||(AniGraph.Group[middleN]==1));
        assert((AniGraph.Group[OppMiddle]==0)||(AniGraph.Group[OppMiddle]==1));

        //enforce consistency of paths
        if (AniGraph.Group[middleN]==1)std::swap(middleN,OppMiddle);

        //then find the first path
        retrievePath(AniGraph,middleN,CurrPath.nodes);
        assert(CurrPath.nodes.size()>=2);
        assert((int)CurrPath.nodes.back()==middleN);

        //get the second path
        std::vector<size_t> otherPath;
        retrievePath(AniGraph,OppMiddle,otherPath);

        assert((int)otherPath.back()==OppMiddle);
        assert(otherPath.size()>=2);

        for (int i=otherPath.size()-2;i>=1;i--)
        {
            int OppI=AniGraph.OppositeNode(otherPath[i]);
            CurrPath.nodes.push_back(OppI);
        }
        return true;
    }

    static bool FindPatternTo(AnisotropicGraph &AniGraph,
                              const size_t startNode,
                              const size_t endNode,
                              const ScalarType max_deviation,
                              const ScalarType penalty_drift,
                              Path &CurrPath)
    {
        DijkstraParam Param;
        Param.sourceGroups.resize(1);
        Param.sourceGroups[0].push_back(startNode);
        Param.maxSink=std::numeric_limits<size_t>::max();
        Param.target=endNode;

        //Param.target
        computePerNodeDijsktra(AniGraph,Param,max_deviation, penalty_drift);
        retrievePath(AniGraph,endNode,CurrPath.nodes);

        //        std::cout<<"Source "<<startNode<<" Target "<<endNode<<std::endl;
        //        std::cout<<"Test"<<CurrPath.nodes.size()<<std::endl;
        if (CurrPath.nodes.size()<=1)return false;
        CurrPath.length=AniGraph.Distance[endNode];
        return true;
    }

    static void FindPatternToSingularities(AnisotropicGraph &AniGraph,
                                           std::vector<size_t> startNodes,
                                           const ScalarType max_deviation,
                                           const ScalarType penalty_drift,
                                           std::vector<Path> &Paths)
    {
        Paths.clear();

        DijkstraParam Param;
        Param.sourceGroups.resize(1);
        Param.sourceGroups[0]=std::vector<size_t>(startNodes.begin(),startNodes.end());
        Param.maxSink=std::numeric_limits<size_t>::max();
        //Param.target
        computePerNodeDijsktra(AniGraph,Param,max_deviation, penalty_drift);

        for (size_t i = AniGraph.nEdgeNodes + AniGraph.nSingNodes; i < AniGraph.Graph.size(); ++i)
        {
            assert(AniGraph.Graph[i]->nodeType == SinkNodeType);
            if (!AniGraph.Graph[i]->Valid)continue;
            Path newPath;
            newPath.length = AniGraph.Distance[i];
            //std::cout << "Cane..." << std::endl;
            retrievePath(AniGraph,i,newPath.nodes);
            if (newPath.nodes.size()<=1)continue;
            Paths.push_back(newPath);
        }
    }

    /**
     * @brief ExtractLoop finds the loop path which starts from and ends to node.
     * @param node
     * @param max_deviation
     * @param penalty_drift
     * @param path
     * @return
     */
    static bool ExtractLoop0(AnisotropicGraph &AniGraph,
                             const size_t node,
                             const ScalarType max_deviation,
                             const ScalarType penalty_drift,
                             Path &path,
                             bool writeMsg=false) {
        // compute distances
        if (writeMsg)
        {
            std::cout << "Computing Dijkstra ... ";
            std::cout.flush();
        }
        int t0 = clock();
        //std::vector<size_t> seeds(1, node);

        DijkstraParam Param;
        Param.sourceGroups.resize(1);
        Param.sourceGroups[0].push_back(node);
        Param.target=node;
        computePerNodeDijsktra(AniGraph,Param,max_deviation, penalty_drift);

        int t1 = clock();

        if (writeMsg)
        {
            std::cout << "done. Time " << t1 - t0;
            // extract reverse path
            std::cout << ". Querying ... ";
            std::cout.flush();
        }

        //std::list<size_t> nodeList;
        path.nodes.clear();
        size_t curr_node = node;
        do {
            //nodeList.push_back(curr_node);
            path.nodes.push_back(curr_node);
            curr_node = AniGraph.Father[curr_node];
        } while (curr_node != std::numeric_limits<size_t>::max() && curr_node != node);
        // copy the path in output
        //path.assign(nodeList.begin(), nodeList.end());
        //int sizeN=path.nodes.size()-1;

        std::cout << "Loop Size " << path.nodes.size()<<std::endl;
        size_t lastN=path.nodes[1];
        assert(lastN>=0);
        path.length=AniGraph.Distance[lastN];

        t1 = clock();

        if (writeMsg)
            std::cout << "done. Time " << t1 - t0;
        std::reverse(path.nodes.begin(),path.nodes.end());
        return curr_node != std::numeric_limits<size_t>::max();
    }

    static ScalarType LoopSelfTangentialDist(AnisotropicGraph &AniGraph,
                                             const std::vector<Path> &Paths,
                                             const ScalarType max_deviation,
                                             const ScalarType penalty_drift,
                                             int N_dir)
    {
        DijkstraParam Param;
        //Param.sourceGroups.resize(Path.nodes.size()*2);
        for (size_t k=0;k<Paths.size();k++)
        {
            for (size_t i=0;i<Paths[k].nodes.size();i++)
            {
                std::vector<size_t> Tangential;
                AniGraph.TangentialNodes(Paths[k].nodes[i],Tangential,N_dir);
                for (size_t j=0;j<Tangential.size();j++)
                {
                    Param.sourceGroups.push_back(std::vector<size_t>());
                    Param.sourceGroups.back().push_back(Tangential[j]);
                }
            }
        }
        int NodeMiddle=computePerNodeDijsktra(AniGraph,Param,max_deviation, penalty_drift);
        if (NodeMiddle==-1)return std::numeric_limits<ScalarType>::max();
        return AniGraph.Distance[NodeMiddle];
    }

    static void getOrthogonalNodes(AnisotropicGraph &AniGraph,
                                   const std::vector<size_t> &Nodes,
                                   int N_dir,std::vector<size_t> &OrthoNodes)
    {
        OrthoNodes.clear();
        for (size_t k=0;k<Nodes.size();k++)
        {
            std::vector<size_t> OrthoN;
            AniGraph.TangentialNodes(Nodes[k],OrthoN,N_dir);
            OrthoNodes.insert(OrthoNodes.end(),OrthoN.begin(),OrthoN.end());
        }
    }

    static void getNodesPositions(AnisotropicGraph &AniGraph,
                                  const std::vector<size_t> &Nodes,
                                  std::vector<CoordType> &NodePositions)
    {
        NodePositions.clear();
        for (size_t k=0;k<Nodes.size();k++)
        {
            CoordType pos=AniGraph.NodePos(Nodes[k]);
            NodePositions.push_back(pos);
        }
    }

    static void ComputeNodesDistances(AnisotropicGraph &AniGraph,
                                      const std::vector<size_t> &Nodes,
                                      const ScalarType max_deviation,
                                      const ScalarType penalty_drift,
                                      std::vector<ScalarType> &NodeDist,
                                      const ScalarType MaxD=std::numeric_limits<ScalarType>::max())
    {
        DijkstraParam Param;

        Param.sourceGroups.resize(1);
        Param.maxDist=MaxD;
        for (size_t i=0;i<Nodes.size();i++)
        {
            size_t NodeI=Nodes[i];
            assert(NodeI>=0);
            assert(NodeI<AniGraph.NumNodes());
            Param.sourceGroups.back().push_back(NodeI);
        }

        //diffuse the loop distance
        if (NodeDist.size()==0)
        {
            computePerNodeDijsktra(AniGraph,Param,max_deviation, penalty_drift);
            NodeDist=std::vector<ScalarType>(AniGraph.Distance.begin(),AniGraph.Distance.end());
        }
        else
        {
            assert(NodeDist.size()==AniGraph.NumNodes());
            std::vector<size_t> withinRadius;
            computePerNodeDijsktra(AniGraph,Param,max_deviation, penalty_drift,&withinRadius,&NodeDist);
            NodeDist.assign(AniGraph.Distance.begin(),AniGraph.Distance.end());
        }
    }

    //    static void ComputeLoopDistances(AnisotropicGraph &AniGraph,
    //                                     const std::vector<Path> &Paths,
    //                                     const ScalarType max_deviation,
    //                                     const ScalarType penalty_drift,
    //                                     std::vector<ScalarType> &NodeDist,
    //                                     const ScalarType MaxD=std::numeric_limits<ScalarType>::max())
    //    {
    //        //NodeDist.resize(AniGraph.Graph.size());
    //        DijkstraParam Param;

    //        Param.sourceGroups.resize(1);
    //        Param.maxDist=MaxD;
    //        for (size_t i=0;i<Paths.size();i++)
    //        {
    //            for (size_t j=0;j<Paths[i].nodes.size();j++)
    //            {
    //                size_t NodeI=Paths[i].nodes[j];
    //                Param.sourceGroups.back().push_back(NodeI);
    //                assert(AniGraph.Graph[NodeI]->Valid);
    //                if (!AniGraph.IsEdgeNode(NodeI))continue;
    //                int OppNode=AniGraph.OppositeNode(NodeI);
    //                assert(OppNode>=0);
    //                assert(OppNode<(int)AniGraph.NumNodes());
    //                Param.sourceGroups.back().push_back(OppNode);
    //            }
    //        }

    //        //diffuse the loop distance
    //        if (NodeDist.size()==0)
    //        {
    //            computePerNodeDijsktra(AniGraph,Param,max_deviation, penalty_drift);
    //            NodeDist=std::vector<ScalarType>(AniGraph.Distance.begin(),AniGraph.Distance.end());
    //        }
    //        else
    //        {
    //            assert(NodeDist.size()==AniGraph.NumNodes());
    //            std::vector<size_t> withinRadius;
    //            computePerNodeDijsktra(AniGraph,Param,max_deviation, penalty_drift,&withinRadius,&NodeDist);
    //            NodeDist.assign(AniGraph.Distance.begin(),AniGraph.Distance.end());
    //        }
    //    }

    //    static void ComputeDistanceFromLoops(AnisotropicGraph &AniGraph,
    //                                     const std::vector<Path> &Paths,
    //                                     const ScalarType max_deviation,
    //                                     const ScalarType penalty_drift,
    //                                     std::vector<ScalarType> &NodeDist,
    //                                     const ScalarType MaxD=std::numeric_limits<ScalarType>::max())
    //    {
    //        //NodeDist.resize(AniGraph.Graph.size());
    //        DijkstraParam Param;

    //        Param.sourceGroups.resize(1);
    //        Param.maxDist=MaxD;
    //        for (size_t i=0;i<Paths.size();i++)
    //        {
    //            for (size_t j=0;j<Paths[i].nodes.size();j++)
    //            {
    //                size_t NodeI=Paths[i].nodes[j];
    //                std::vector<int> TNodes;
    //                LoopFunctions<MeshType>::GetTangentNodes(AniGraph,NodeI,TNodes);
    //                for (size_t )
    //                Param.sourceGroups.back().push_back(NodeI);
    //                assert(AniGraph.Graph[NodeI]->Valid);
    //                if (!AniGraph.IsEdgeNode(NodeI))continue;
    //                int OppNode=AniGraph.OppositeNode(NodeI);
    //                assert(OppNode>=0);
    //                assert(OppNode<(int)AniGraph.NumNodes());
    //                Param.sourceGroups.back().push_back(OppNode);
    //            }
    //        }

    //        //diffuse the loop distance
    //        if (NodeDist.size()==0)
    //        {
    //            computePerNodeDijsktra(AniGraph,Param,max_deviation, penalty_drift);
    //            NodeDist=std::vector<ScalarType>(AniGraph.Distance.begin(),AniGraph.Distance.end());
    //        }
    //        else
    //        {
    //            assert(NodeDist.size()==AniGraph.NumNodes());
    //            std::vector<size_t> withinRadius;
    //            computePerNodeDijsktra(AniGraph,Param,max_deviation, penalty_drift,&withinRadius,&NodeDist);
    //            NodeDist.assign(AniGraph.Distance.begin(),AniGraph.Distance.end());
    //        }
    //    }

    static void ComputeDistanceFromNodes(AnisotropicGraph &AniGraph,
                                         const std::vector<int> &Nodes,
                                         const ScalarType max_deviation,
                                         const ScalarType penalty_drift,
                                         std::vector<ScalarType> &NodeDist,
                                         const ScalarType MaxD=std::numeric_limits<ScalarType>::max())
    {
        //NodeDist.resize(AniGraph.Graph.size());
        DijkstraParam Param;

        Param.sourceGroups.resize(1);
        Param.maxDist=MaxD;
        for (size_t i=0;i<Nodes.size();i++)
        {
            int NodeI=Nodes[i];
            Param.sourceGroups.back().push_back(NodeI);
            assert(AniGraph.Graph[NodeI]->Valid);
        }

        //diffuse the loop distance
        if (NodeDist.size()==0)
        {
            computePerNodeDijsktra(AniGraph,Param,max_deviation, penalty_drift);
            NodeDist=std::vector<ScalarType>(AniGraph.Distance.begin(),AniGraph.Distance.end());
        }
        else
        {
            assert(NodeDist.size()==AniGraph.NumNodes());
            std::vector<size_t> withinRadius;
            computePerNodeDijsktra(AniGraph,Param,max_deviation, penalty_drift,&withinRadius,&NodeDist);
            NodeDist.assign(AniGraph.Distance.begin(),AniGraph.Distance.end());
        }
    }

    //    static void ComputeSingleLoopDistances(AnisotropicGraph &AniGraph,
    //                                           const std::vector<Path> &Paths,
    //                                           const ScalarType max_deviation,
    //                                           const ScalarType penalty_drift,
    //                                           const Path &NewLoop,
    //                                           std::vector<ScalarType> &NodeDist)
    //    {
    //        NodeDist.resize(NewLoop.nodes.size());
    //        DijkstraParam Param;
    //        //set the other nodes as
    //        Param.sourceGroups.resize(1);
    //        for (size_t i=0;i<Paths.size();i++)
    //        {
    //            for (size_t j=0;j<Paths[i].nodes.size();j++)
    //            {
    //                size_t NodeI=Paths[i].nodes[j];
    //                Param.sourceGroups.back().push_back(NodeI);
    //                int OppNode=AniGraph.OppositeNode(NodeI);
    //                assert(OppNode>=0);
    //                assert(OppNode<AniGraph.NumNodes());
    //                Param.sourceGroups.back().push_back(OppNode);
    //            }
    //        }

    //        //diffuse the loop distance
    //        computePerNodeDijsktra(AniGraph,Param,max_deviation, penalty_drift);
    //        //then return per-node distance
    //        NodeDist=std::vector<ScalarType>(AniGraph.Distance.begin(),AniGraph.Distance.end());
    //    }


    static bool ExtractLoop(AnisotropicGraph &AniGraph,
                            const size_t node,
                            const ScalarType max_deviation,
                            const ScalarType penalty_drift,
                            Path &path,
                            bool writeMsg=false) {
        // compute distances
        if (writeMsg)
        {
            std::cout << "Computing Dijkstra ... ";
            std::cout.flush();
        }
        int t0 = clock();
        //std::vector<size_t> seeds(1, node);

        DijkstraParam Param;
        Param.sourceGroups.resize(2);
        int OppN=AniGraph.OppositeNode(node);
        assert(OppN!=-1);
        assert(AniGraph.Graph[node]->Valid);
        assert(AniGraph.Graph[OppN]->Valid);
        Param.sourceGroups[0].push_back(node);
        Param.sourceGroups[1].push_back(OppN);

        //Param.target=node;
        int middleN=computePerNodeDijsktra(AniGraph,Param,max_deviation, penalty_drift);

        if (middleN==-1)return false;

        int t1 = clock();
        //std::cout<<"Tracing Time "<<(t1-t0)/(ScalarType)CLOCKS_PER_SEC<<std::endl;

        if (writeMsg)
        {
            std::cout << "done. Time (msec) :" << (t1 - t0)/( CLOCKS_PER_SEC / 1000 )<<std::endl;
            // extract reverse path
            std::cout << "Retrieving path ... ";
            std::cout.flush();
        }

        path.nodes.clear();

        //get the opposite node
        int OppMiddle=AniGraph.OppositeNode(middleN);
        assert(OppMiddle!=-1);
        //set the lenght
        path.length=AniGraph.Distance[middleN];
        path.length+=AniGraph.Distance[OppMiddle];

        //check group consistency
        assert(AniGraph.Group[middleN]!=AniGraph.Group[OppMiddle]);
        assert((AniGraph.Group[middleN]==0)||(AniGraph.Group[middleN]==1));
        assert((AniGraph.Group[OppMiddle]==0)||(AniGraph.Group[OppMiddle]==1));

        //enforce consistency of paths
        if (AniGraph.Group[middleN]==1)std::swap(middleN,OppMiddle);

        //then find the first path
        std::cout.flush();
        retrievePath(AniGraph,node,middleN,path.nodes);
        //std::cout << "size ... " << path.nodes.size() << std::endl;
        assert(path.nodes.size()>=2);
        assert((int)path.nodes.back()==middleN);
        assert(path.nodes[0]==node);

        //get the second path
        std::vector<size_t> otherPath;
        retrievePath(AniGraph,OppN,OppMiddle,otherPath);
        //std::cout << "size ... " << otherPath.size() << std::endl;
        //std::cout.flush();

        assert((int)otherPath.back()==OppMiddle);
        assert((int)otherPath[0]==OppN);
        assert(otherPath.size()>=2);
        for (int i=otherPath.size()-2;i>=1;i--)
        {
            int OppI=AniGraph.OppositeNode(otherPath[i]);
            path.nodes.push_back(OppI);
        }
        // copy the path in output
        int t2 = clock();

        if (writeMsg)
            std::cout << "done. Time Time (msec): " << (t2 - t1)/( CLOCKS_PER_SEC / 1000 )<<std::endl;

        //std::reverse(path.nodes.begin(),path.nodes.end());
        return true;//curr_node != std::numeric_limits<size_t>::max();
    }


    static void GetTangentNodes(AnisotropicGraph &AniGraph,
                                const int IndexNode,
                                int &OrthoSeed0,
                                int &OrthoSeed1)
    {

        std::vector<size_t> Tangential;
        AniGraph.TangentialNodes(IndexNode,Tangential);
        OrthoSeed0=-1;
        OrthoSeed1=-1;
        if (Tangential.size()==2)
        {
            if ((AniGraph.IsValid(Tangential[0]))&&
                    (AniGraph.IsValid(Tangential[1])))
            {
                OrthoSeed0=Tangential[0];
                OrthoSeed1=Tangential[1];
                return;
            }
        }
        int OppNode=AniGraph.OppositeNode(IndexNode);
        AniGraph.TangentialNodes(OppNode,Tangential);
        if (Tangential.size()==2)
        {

            if ((AniGraph.IsValid(Tangential[0]))&&
                    (AniGraph.IsValid(Tangential[1])))
            {
                OrthoSeed0=Tangential[0];
                OrthoSeed1=Tangential[1];
            }
        }
    }

    static void MoveOneStep(AnisotropicGraph &AniGraph,
                            const int IndexNode,
                            ScalarType max_deviation,
                            ScalarType penalty_drift,
                            int &IndexSucc)
    {
        IndexSucc=-1;
        ScalarType bestDist=std::numeric_limits<ScalarType>::max();
        //then move of one step in tangential dir
        for(NeighborsIterator it = AniGraph.Graph[IndexNode]->neighbors.begin();
            it != AniGraph.Graph[IndexNode]->neighbors.end(); ++it)
        {
            //discard if it is not valid
            if (!it->Valid)continue;

            int IndexNeigh=it->nodeID;
            assert(IndexNeigh < AniGraph.Graph.size());

            if (!AniGraph.IsValid(IndexNeigh))continue;

            // discard neigbors outside the cone
            if (it->angle > max_deviation)
                continue;

            // compute the new distance to the current neighbor
            ScalarType nextDistance = it->geoDistance(max_deviation, penalty_drift);
            if (nextDistance>bestDist)continue;
            bestDist=nextDistance;

            IndexSucc=IndexNeigh;
        }
    }

    static void DistanceStats(AnisotropicGraph &AniGraph,
                              const ScalarType max_deviation,
                              const ScalarType penalty_drift,
                              ScalarType &MinDist,
                              ScalarType &MaxDist,
                              ScalarType &AvDist)
    {
        MinDist=std::numeric_limits<ScalarType>::max();
        MaxDist=0;
        AvDist=0;
        size_t NumAv=0;
        for (size_t i=0;i<AniGraph.NumNodes();i++)
        {
            for(NeighborsIterator it = AniGraph.Graph[i]->neighbors.begin();
                it != AniGraph.Graph[i]->neighbors.end(); ++it)
            {
                // compute the new distance to the current neighbor
                ScalarType testDistance = it->geoDistance(max_deviation, penalty_drift);
                MinDist=std::min(MinDist,testDistance);
                MaxDist=std::min(MaxDist,testDistance);
                AvDist+=testDistance;
                NumAv++;
            }
        }
        AvDist/=NumAv;
    }
    //    static void GetNighNodes(AnisotropicGraph &AniGraph,
    //                             const int &IndexNode,
    //                             const int &NeighStep,
    //                             std::vector<size_t> &NeighNodes)
    //    {
    //        assert(NeighStep>0);
    //        NeighNodes.clear();
    //        NeighNodes.push_back(IndexNode);

    //        //for each step
    //        int PrevStep=0;
    //        for (size_t s=0;s<NeighStep;s++)
    //        {
    //            std::vector<size_t> NextStepNodes;
    //            for (size_t n=PrevStep;n<NeighNodes.size();n++)
    //            {
    //                size_t currN=NeighNodes[n];
    //                for(NeighborsIterator it = AniGraph.Graph[currN]->neighbors.begin();
    //                    it != AniGraph.Graph[currN]->neighbors.end(); ++it)
    //                {
    //                    int IndexNeigh=it->nodeID;
    //                    NextStepNodes.push_back(IndexNeigh);
    //                }
    //            }
    //            std::sort(NextStepNodes.begin(), NextStepNodes.end());
    //            std::vector<size_t>::iterator last = std::unique(NextStepNodes.begin(), NextStepNodes.end());
    //            NextStepNodes.erase(last, NextStepNodes.end());
    //            //then set the new starting poin
    //            PrevStep=NeighNodes.size();

    //            //and add new nodes
    //            NeighNodes.insert(NeighNodes.end(),NextStepNodes.begin(),NextStepNodes.end());
    //        }
    //    }
};
#endif //ANISOTROPIC_GEODESIC_QUERY
