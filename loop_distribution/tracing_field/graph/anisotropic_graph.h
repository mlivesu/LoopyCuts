#ifndef ANISOTROPIC_GEODESIC_GRAPH
#define ANISOTROPIC_GEODESIC_GRAPH

#include "base_geodesic_node.h"
#include "edge_geodesic_node.h"
#include "singulatity_geodesic_node.h"
#include "sink_geodesic_node.h"
#include "neigh_info.h"
#include "path_geodesic_node.h"
#include <tracing_field/stats_collector.h>
#include <vcg/complex/algorithms/update/color.h>
#ifndef NO_TRACING_OPENGL
#include <wrap/gl/addons.h>
#endif
#include <vcg/complex/algorithms/parametrization/tangent_field_operators.h>

/**
 * @brief The AnisotropicGraph class
 */
template < class MeshType >
class AnisotropicGraph {
    typedef typename MeshType::CoordType    CoordType;
    typedef typename MeshType::VertexType   VertexType;
    typedef typename MeshType::FaceType     FaceType;
    typedef typename MeshType::ScalarType   ScalarType;

    typedef EdgeGeodesicNode<CoordType> EdgeGeodesicNode;
    typedef SingularityGeodesicNode<CoordType> SingularityGeodesicNode;
    typedef SinkGeodesicNode<CoordType> SinkGeodesicNode;

    typedef NeighborInfo<ScalarType> NeighborInfo;
    typedef BaseGeodesicNode<CoordType> BaseGeodesicNode;

public:
    typedef Path<ScalarType> Path;
    typedef typename BaseGeodesicNode::NeighborsIterator NeighborsIterator;

    /**
     * @brief AnisotropicGraph Constructor.
     * @param _mesh
     */
    AnisotropicGraph(MeshType &_mesh,
                     StatsCollector<ScalarType> &_Stats) : nEdgeNodes(0),
        nSingNodes(0),
        mesh(_mesh),
        Stats(_Stats),
        worst_connection(std::numeric_limits<ScalarType>::max()),
        best_connection(std::numeric_limits<ScalarType>::max()),
        max_angle(45)
    { }

    /**
     * @brief ~AnisotropicGraph Destructor.
     */
    ~AnisotropicGraph() {
        clearGraph();
    }


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
    
    typedef std::list<size_t>                               NodeList;
    typedef typename NodeList::iterator                     NodeListIterator;
    typedef typename NodeList::const_iterator               NodeListConstIterator;
    size_t                                                  nEdgeNodes;         ///< Number of edge nodes.
    size_t                                                  nSingNodes;         ///< Number of singularity nodes.
    std::deque<BaseGeodesicNode*>                           Graph;              ///< The actual graph.

private:
    
    // private data members *******************************************************************************************************************************************
    MeshType                                               &mesh;               ///< The input mesh.
    StatsCollector<ScalarType>                             &Stats;              ///< Class to collect the statistics
    ScalarType                                              worst_connection;
    ScalarType                                              best_connection;
    ScalarType                                              max_angle;          ///< Maximum deviation allowed.
    std::vector<std::vector<std::vector<NodeList> > >       PerEdgeNodes;       ///< Per-face, per-edge and per-direction nodes (at each sample point).
    std::map<CoordType, NodeList>                           CoordNodes;         ///< Map of the nodes that share the same coordinates.
    std::vector<NodeList>                                   PerVertNodes;       ///< Per-vertex nodes (singularities and sinks).

    // private methods ************************************************************************************************************************************************

    /**
     * @brief clearGraph deletes the graph.
     */
    void clearGraph() {
        for (size_t i = 0; i < Graph.size(); ++i)
            delete Graph[i];
        Graph.clear();
        worst_connection = std::numeric_limits<ScalarType>::max();
        best_connection = std::numeric_limits<ScalarType>::max();
        nEdgeNodes = 0;
        nSingNodes = 0;
        PerEdgeNodes.clear();
        CoordNodes.clear();
        PerVertNodes.clear();
    }

    /**
     * @brief edgeSize computes the minimum, maximum and average edge lengths.
     * @param minSize
     * @param maxSize
     * @param avgSize
     */
    void edgeSize(ScalarType &minSize,
                  ScalarType &maxSize,
                  ScalarType &avgSize) {
        avgSize = maxSize = ScalarType(0);
        minSize = std::numeric_limits<ScalarType>::max();
        size_t numE = 0;
        for (size_t i = 0; i < mesh.face.size(); ++i)
            for (int j = 0; j < mesh.face[i].VN(); ++j)
                if (vcg::face::IsBorder(mesh.face[i], j) || mesh.face[i].V0(j) < mesh.face[i].V1(j)) {
                    ScalarType curr_size = (mesh.face[i].P1(j) - mesh.face[i].P0(j)).Norm();
                    if (curr_size < minSize)
                        minSize = curr_size;
                    if (curr_size > maxSize)
                        maxSize = curr_size;
                    avgSize += curr_size;
                    ++numE;
                }
        avgSize /= static_cast<ScalarType>(numE);
    }

    /**
     * @brief samplePos computes a new point on an edge by interpolation.
     * @param f
     * @param edgeIndex
     * @param subDivision
     * @param curr_pos
     * @return
     */
    CoordType samplePos(const FaceType &f,
                        const int edgeIndex,
                        const size_t subDivision,
                        const size_t curr_pos) {
        assert(edgeIndex >= 0 && edgeIndex < f.VN());
        assert(subDivision > 0);
        assert(curr_pos <= subDivision - 1);
        // interpolate in the current pos
        ScalarType alpha = curr_pos / static_cast<ScalarType>(subDivision);
        return f.cP0(edgeIndex) * (ScalarType(1) - alpha) + f.cP1(edgeIndex) * alpha;
    }

    /**
     * @brief addEdgeSample
     * @param indexF
     * @param indexE
     * @param pos
     */
    //THIS MUST BE EXTENDED FOR ARBITRARY DIR NUMBERS
    void addEdgeSample(const size_t indexF,
                       const int indexE,
                       const CoordType &pos,
                       bool IsBorderE)
    {
        std::vector<size_t> EdgeNodesIndex;

        //add all 4 directions of the cross field
        for (int dir = 0; dir < 4; ++dir)
        {
            Graph.push_back(new EdgeGeodesicNode(indexF, indexE, dir, pos,IsBorderE));

            EdgeNodesIndex.push_back(Graph.size()-1);
            //EdgeNodesPtr.push_back((EdgeGeodesicNode*) Graph.back());

            PerEdgeNodes[indexF][indexE][dir].push_back(Graph.size() - 1);
        }

        //then initialize opposite of the same face
        for (size_t i=0;i<EdgeNodesIndex.size();i++)
        {
            int IndexN=EdgeNodesIndex[i];

            assert(Graph[IndexN]->nodeType == EdgeNodeType);
            EdgeGeodesicNode *edgeN = static_cast<EdgeGeodesicNode*>(Graph[IndexN]);

            edgeN->OppositeSameFID=EdgeNodesIndex[(i+2)%4];
        }
        nEdgeNodes += 4;
    }

    /**
     * @brief sampleMesh computes edge-nodes samples.
     * @param subDivisionStep
     */
    void sampleMesh(const int NumSubdivisions,
                    bool UseQ=true)
    {
        clearGraph();
        //update quality per face if required
        if (UseQ)
        {
            for (size_t i=0;i<mesh.face.size();i++)
            {
                //update the quality
                mesh.face[i].Q()=vcg::QualityRadii(mesh.face[i].P(0),mesh.face[i].P(1),mesh.face[i].P(2));
                //mesh.face[i].Q()=pow(mesh.face[i].Q(),2);
                //compute the multiplier
                mesh.face[i].Q()=1.0+((1.0-mesh.face[i].Q())*(NumSubdivisions*5));
            }
        }

        // allocate space for references
        PerEdgeNodes.resize(mesh.face.size());
        for (size_t i = 0; i < mesh.face.size(); ++i)
            PerEdgeNodes[i].resize(mesh.face[i].VN(), std::vector<NodeList>(4));

        size_t Avg=0;
        size_t SamplStep=0;
        // initialize per edge samples
        for (size_t i = 0; i < mesh.face.size(); ++i) {
            for (int j = 0; j < mesh.face[i].VN(); ++j)
            {
                bool IsBorderE=(mesh.face[i].FFp(j)==&mesh.face[i]);

                // edges must be sampled only once check the order and the fact that is not a border edge
                if ((vcg::tri::Index(mesh, mesh.face[i].V0(j)) >= vcg::tri::Index(mesh, mesh.face[i].V1(j)))&&
                        (!IsBorderE))
                    continue;

                // compute samples
                size_t multiplier=1;
                if (UseQ)
                {
                    ScalarType Q0=mesh.face[i].Q();
                    ScalarType Q1=mesh.face[i].FFp(j)->Q();
                    multiplier*=std::max(Q0,Q1);
                }
                //ScalarType edgeSize = (mesh.face[i].V1(j)->P() - mesh.face[i].V0(j)->P()).Norm();
                //size_t numSub = (2 + floor(edgeSize / subDivisionStep + 0.5)) * 2;
                int numSub = multiplier*2;//(floor(edgeSize / subDivisionStep + 0.5)) * 2;

                if (numSub<NumSubdivisions*2)
                    numSub=NumSubdivisions*2;

                SamplStep++;
                for (int k = 1; k < numSub; k += 2) {
                    Avg++;
                    CoordType pos = samplePos(mesh.face[i], j, numSub, k);
                    addEdgeSample(i, j, pos,IsBorderE);

                    //if is not border then add the node on opposite face
                    if (IsBorderE)continue;

                    addEdgeSample(vcg::tri::Index(mesh, mesh.face[i].FFp(j)), mesh.face[i].FFi(j), pos,IsBorderE);
                }
            }
        }
        std::cout << "Sampled in Average " << (float)Avg/(float)SamplStep << " samples per Edge " << std::endl;
    }

    /**
     * @brief getPosSurroundingE returns the surrounding edges of the given edge outside the given face.
     * @param indexF
     * @param indexE
     * @param surroundEdges
     */
    void getPosSurroundingE(const size_t indexF,
                            const int indexE,
                            std::list<vcg::face::Pos<FaceType> > &surroundEdges)
    {
        surroundEdges.clear();
        int indexE0=(indexE+1)%3;
        int indexE1=(indexE+2)%3;

        FaceType *F0=mesh.face[indexF].FFp(indexE0);
        FaceType *F1=mesh.face[indexF].FFp(indexE1);
        int Eopp0=mesh.face[indexF].FFi(indexE0);
        int Eopp1=mesh.face[indexF].FFi(indexE1);
        //check if it is not a border edge
        if (Eopp0!=-1) surroundEdges.push_back(vcg::face::Pos<FaceType>(F0,Eopp0));
        if (Eopp1!=-1) surroundEdges.push_back(vcg::face::Pos<FaceType>(F1,Eopp1));
    }

    /**
     * @brief getPossibleConnections returns the nodes that can be connected to the given one.
     * @param indexNode
     * @param nodes
     */
    void getPossibleConnections(const size_t indexNode,
                                std::list<size_t> &nodes) {
        //check if the node is a vertex node
        assert(Graph[indexNode]->nodeType == EdgeNodeType);
        EdgeGeodesicNode *edgeN = static_cast<EdgeGeodesicNode*>(Graph[indexNode]);

        // get the other edges
        std::list<vcg::face::Pos<FaceType> > surroundEdges;
        getPosSurroundingE(edgeN->faceID, edgeN->edgeIdx, surroundEdges);

        // then given the direction of the current face
        // get the direction on the closest face
        typename std::list<vcg::face::Pos<FaceType> >::iterator it;
        for (it = surroundEdges.begin(); it != surroundEdges.end(); ++it) {
            // get face and edge on other side
            size_t indexF1 = vcg::tri::Index(mesh, it->F());
            int e1 = it->E();
            // get the direction it has to follow
            int oppDir = vcg::tri::CrossField<MeshType>::FollowDirection(mesh.face[edgeN->faceID], *(it->F()), edgeN->M4Dir);
            // collect the possible connections
            nodes.insert(nodes.end(), PerEdgeNodes[indexF1][e1][oppDir].begin(), PerEdgeNodes[indexF1][e1][oppDir].end());
        }
    }

    /**
     * @brief isCompatible checks if two nodes are compatible and computes the weight.
     * @param node0
     * @param node1
     * @param currAngle
     * @param length
     * @return
     */
    bool isCompatible(const size_t node0,
                      const size_t node1,
                      ScalarType &currAngle,
                      ScalarType &length) {
        assert(node0 < Graph.size());
        assert(node1 < Graph.size());

        assert(Graph[node0]->nodeType == EdgeNodeType);
        EdgeGeodesicNode *edgeN0 = static_cast<EdgeGeodesicNode*>(Graph[node0]);
        assert(Graph[node1]->nodeType == EdgeNodeType);
        EdgeGeodesicNode *edgeN1 = static_cast<EdgeGeodesicNode*>(Graph[node1]);

        //START CHECK DIRECTION
        assert(edgeN0->M4Dir >= 0 && edgeN0->M4Dir < 4);
        assert(edgeN1->M4Dir >= 0 && edgeN1->M4Dir < 4);
        assert((edgeN0->faceID) < mesh.face.size());
        assert((edgeN1->faceID) < mesh.face.size());
        assert(vcg::tri::CrossField<MeshType>::FollowDirection(mesh.face[edgeN0->faceID], mesh.face[edgeN1->faceID], edgeN0->M4Dir) == edgeN1->M4Dir);
        //END CHECK DIRECTION

        // get the direction of the current face
        CoordType dir1 = edgeN1->pos - edgeN0->pos;
        // compute its length
        length = dir1.Norm();

        // compute the angle
        dir1.Normalize();
        currAngle = fabs(vcg::Angle(getDirection(node0), dir1) * 180.0 / M_PI);

        // check if the angle is within the interval
        if (currAngle > max_angle)
            return false;
        //    if (fabs(vcg::Angle(dir1, getDirection(node1)) * 180.0 / M_PI) > max_angle)
        //        return false;
        return true;
    }

    /**
     * @brief getConnections returns the nodes that can be connected to the given one, together with the corresponding weights.
     * @param indexNode
     * @param neighbors
     */
    void getConnections(const size_t indexNode,
                        std::list<NeighborInfo > &neighbors) {
        neighbors.clear();

        // get the neighboring possible nodes
        std::list<size_t> nodes;
        getPossibleConnections(indexNode, nodes);

        // filter out the ones that deviates too much from chosen direction
        ScalarType angle, length;
        for (std::list<size_t>::iterator it = nodes.begin(); it != nodes.end(); ++it)
        {
            bool compatible=isCompatible(indexNode, *it, angle, length);
            //if (isCompatible(indexNode, *it, angle, length))
            if (compatible)
                neighbors.push_back(NeighborInfo(*it, angle, length));
        }
    }

    /**
     * @brief initConnections initializes edge-nodes' connections.
     */
    void initConnections() {
        for (size_t i = 0; i < nEdgeNodes; ++i) {
            assert(Graph[i]->nodeType == EdgeNodeType);
            EdgeGeodesicNode *edgeN = static_cast<EdgeGeodesicNode*>(Graph[i]);
            getConnections(i, edgeN->neighbors);
        }
    }

    /**
     * @brief initCoordTable
     */
    void initCoordTable() {
        CoordNodes.clear();
        for (size_t i = 0; i < Graph.size(); ++i)
            CoordNodes[Graph[i]->pos].push_back(i);
    }

    int ComputeOppositeEdgeNode(size_t NodeI,int NDir=2)
    {
        assert(IsEdgeNode(NodeI));

        EdgeGeodesicNode *edgeN = static_cast<EdgeGeodesicNode*>(Graph[NodeI]);

        //get the face, edge and M4Dir
        int FIndex=edgeN->faceID;
        int EIndex=edgeN->edgeIdx;
        size_t M4Dir=edgeN->M4Dir;

        //get opposite direction in current face
        size_t OppDir0=(M4Dir+NDir)%(NDir*2);

        //get the face on the opposite side
        FaceType *F0=&mesh.face[FIndex];
        FaceType *Fopp=F0->FFp(EIndex);
        int EOppIndex=F0->FFi(EIndex);

        size_t FOppIndex=vcg::tri::Index(mesh,Fopp);
        assert((int)FOppIndex!=FIndex);//not on the border

        //get opposite direction in adjacent face
        size_t M4OppDir=vcg::tri::CrossField<MeshType>::FollowDirection(*F0,*Fopp,OppDir0);

        //then get the ones on the opposite edge that have the same pos
        CoordType testPos0=edgeN->pos;

        //then check the one with the same position
        for (NodeListConstIterator NListIte=PerEdgeNodes[FOppIndex][EOppIndex][M4OppDir].begin();
             NListIte!=PerEdgeNodes[FOppIndex][EOppIndex][M4OppDir].end();
             NListIte++)
        {
            size_t CurrIndex=(*NListIte);
            assert(IsEdgeNode(CurrIndex));
            EdgeGeodesicNode *edgeN1 = static_cast<EdgeGeodesicNode*>(Graph[CurrIndex]);
            CoordType testPos1=edgeN1->pos;
            if (testPos0==testPos1)
            {
                return(CurrIndex);
            }
        }
        return -1;
    }

    void InitOppositeNodes(int NDir=2)
    {
        NoOppNodes.clear();
        for (size_t i = 0; i < Graph.size(); ++i)
        {
            if (!IsEdgeNode(i))continue;
            if (IsBorderNode(i))continue;
            int OppNode=ComputeOppositeEdgeNode(i,NDir);
            if (OppNode==-1)
            {
                NoOppNodes.push_back(i);
                std::cout << "WARNING NO OPPOSITE NODE " << std::endl;
            }
            EdgeGeodesicNode *edgeN = static_cast<EdgeGeodesicNode*>(Graph[i]);
            edgeN->OppositeID=OppNode;
        }
    }

public:

    size_t expectedValence(const int &vIndex)
    {return vcg::tri::CrossField<MeshType>::expectedValence(mesh,mesh.vert[vIndex]);}
    //{return(mesh.expectedValence(mesh.vert[vIndex]));}
    //{return(expectedValence(mesh.vert[vIndex]));}

    bool IsSingularVert(const int &vIndex)const
    {return vcg::tri::CrossField<MeshType>::IsSingular(mesh,mesh.vert[vIndex]);}
    //{return(mesh.IsSingularVert(mesh.vert[vIndex]));}

    //{return(IsSingularVert(mesh.vert[vIndex]));}

private:

    ///**
    //     * @brief findSeparatrices
    //     * @param vPos
    //     * @param directions
    //     */

    std::vector<vcg::Triangle3<ScalarType> > WrongTris;
    std::vector<size_t> WrongNodes;
    std::vector<size_t> NoOppNodes;


    void ParametrizeStar(vcg::face::Pos<FaceType> CentralPos)
    {
        //get the element in the star
        std::vector<VertexType*> vertexVec;
        vcg::face::VVOrderedStarFF(CentralPos,vertexVec);

        //find the sum of angles
        ScalarType intervalAngle=2*M_PI/(ScalarType)vertexVec.size();

        VertexType *CenterV=CentralPos.V();

        CenterV->T().P().X()=0;
        CenterV->T().P().Y()=0;

        ScalarType CurrAngle=0;
        for (size_t i = 0; i < vertexVec.size(); ++i)
        {
            vertexVec[i]->T().P().X()=cos(CurrAngle);
            vertexVec[i]->T().P().Y()=sin(CurrAngle);
            CurrAngle+=intervalAngle;
        }
    }


    /**
     * @brief testSeparatrices
     */
    void testSeparatrices() {
        WrongTris.clear();
        size_t valence, vIndex;
        vcg::face::Pos<FaceType> currP;
        //std::list<CoordType> directions,norms;

        vcg::tri::UpdateFlags<MeshType>::VertexClearS(mesh);
        for (size_t i = 0; i < mesh.face.size(); ++i) {
            if (mesh.face[i].IsD())
                continue;
            for (int j = 0; j < mesh.face[i].VN(); ++j)
            {
                if (mesh.face[i].V(j)->IsB())
                    continue;
                if (mesh.face[i].V(j)->IsS())
                    continue;
                vIndex = vcg::tri::Index(mesh,mesh.face[i].V(j));
                mesh.vert[vIndex].SetS();
                valence = mesh.expectedValence(mesh.vert[vIndex]);
                currP.Set(&mesh.face[i], j, mesh.face[i].V(j));
                //bool isOK=findSeparatrices(currP, directions);
                std::vector<FaceType*> Faces;
                std::vector<CoordType> directions;
                size_t sampledDir=vcg::tri::CrossField<MeshType>::FindSeparatrices(currP,directions,Faces,WrongTris,valence);
                //assert(valence == directions.size());
                if (sampledDir!=directions.size()){
                    std::cout << "Sampled " << sampledDir << "directions instead of " << valence << " CORRECTED!" <<std::endl;
                }
            }
        }
    }

    /**
     * @brief addSingularityNodes creates a new node for each direction and adds it to the graph.
     * @param vIndex
     * @param directions
     * @param newNodes
     */
    void addSingularityNodes(const size_t vIndex,
                             const std::vector<CoordType> &directions,
                             const std::vector<FaceType*> &Faces,
                             std::vector<size_t> &newNodes)
    {
        assert(directions.size()==Faces.size());
        newNodes.clear();
        for (size_t i=0;i<directions.size();i++)
        {
            size_t IndexF=vcg::tri::Index(mesh,Faces[i]);
            // add to the graph
            Graph.push_back(new SingularityGeodesicNode(vIndex, mesh.vert[vIndex].P(), directions[i],IndexF));
            newNodes.push_back(Graph.size() - 1);
            // update auxiliary structures
            PerVertNodes[vIndex].push_back(Graph.size() - 1);
        }
        nSingNodes += directions.size();
    }

    /**
     * @brief faceEdgeNodes returns the nodes on a given edge belonging to a given face.
     * @param indexFace
     * @param indexEdge
     * @param nodes
     */
    void faceEdgeNodes(const size_t indexFace,
                       const int indexEdge,
                       std::map<CoordType, NodeList> &nodes) {
        for (size_t i = 0; i < PerEdgeNodes[indexFace][indexEdge].size(); ++i)
            for (NodeListIterator it = PerEdgeNodes[indexFace][indexEdge][i].begin(); it != PerEdgeNodes[indexFace][indexEdge][i].end(); ++it)
                nodes[Graph[*it]->pos].push_back(*it);
    }

public:

    void faceEdgeNodesDir(const size_t indexFace,
                          const int indexEdge,
                          const int M4Dir,
                          std::vector<size_t> &nodesI)const
    {
        //first get the closest face direction
        for (size_t i = 0; i < PerEdgeNodes[indexFace][indexEdge].size(); ++i)
            for (NodeListConstIterator it = PerEdgeNodes[indexFace][indexEdge][i].begin();
                 it != PerEdgeNodes[indexFace][indexEdge][i].end(); ++it)
            {
                size_t IndexN=(*it);
                EdgeGeodesicNode *edgeN = static_cast<EdgeGeodesicNode*>(Graph[IndexN]);
                if (edgeN->M4Dir!=M4Dir)continue;
                nodesI.push_back(IndexN);
            }
    }


    void faceEdgeNodesDirCoord(const size_t indexFace,
                               const int indexEdge,
                               const CoordType Dir,
                               std::vector<size_t> &nodesI)const
    {
        ScalarType MaxDot=-1;
        int BestDir=-1;
        for (size_t i=0;i<4;i++)
        {
            CoordType TestDir=vcg::tri::CrossField<MeshType>::CrossVector(mesh.face[indexFace], i);
            if ((TestDir*Dir)<MaxDot)continue;
            MaxDot=TestDir*Dir;
            BestDir=i;
        }
        assert(BestDir!=-1);
        faceEdgeNodesDir(indexFace,indexEdge,BestDir,nodesI);
    }

    void faceNodesDir(const size_t indexFace,
                      const int M4Dir,
                      std::vector<size_t> &nodesI)
    {
        faceEdgeNodesDir(indexFace,0,M4Dir,nodesI);
        faceEdgeNodesDir(indexFace,1,M4Dir,nodesI);
        faceEdgeNodesDir(indexFace,2,M4Dir,nodesI);
    }

    //    void ReduceWeightLinks(std::vector<size_t> &IndexNodes,
    //                           ScalarType InsideW=0.00000001)
    //    {
    //        std::set<size_t> NodeSet(IndexNodes.begin(),IndexNodes.end());
    //         for (size_t i=0;i<Graph.size();i++)
    //        {
    //             for (NeighborsIterator NeighIte=Graph[i]->neighbors.begin();
    //                 NeighIte!=Graph[i]->neighbors.end();NeighIte++)
    //            {
    //                //if (NodeSet.count(i)==0)continue;
    //                if (NodeSet.count(NeighIte->nodeID)==0)continue;

    //                (*NeighIte).Weight=InsideW;
    //            }
    //        }
    //    }
private:



    /**
     * @brief faceEdgeNodes returns the nodes on the edges given in posVec.
     * @param posVec
     * @param nodes
     */
    void faceEdgeNodes(const std::list<vcg::face::Pos<FaceType> > &posVec,
                       std::map<CoordType, NodeList> &nodes) {
        nodes.clear();
        for (typename std::list<vcg::face::Pos<FaceType> >::const_iterator it = posVec.begin(); it != posVec.end(); ++it)
            faceEdgeNodes(vcg::tri::Index(mesh, it->F()), it->E(), nodes);
    }



    /**
     * @brief getPosSurroundingESing returns the outside edges in the one-ring of the given singular vertex.
     * @param indexF
     * @param indexE
     * @param surroundEdges
     */
    void getPosSurroundingESing(const vcg::face::Pos<FaceType> &singPos,
                                std::list<vcg::face::Pos<FaceType> > &surroundEdges)
    {
        //    assert(indexF < mesh.face.size());
        //    assert(indexE >= 0 && indexE < mesh.face[indexF].VN());

        assert(singPos.F()->V(singPos.E())==singPos.V());
        int indexF=vcg::tri::Index(mesh,singPos.F());
        int indexE=singPos.E();

        surroundEdges.clear();

        vcg::face::Pos<FaceType> start_pos(&mesh.face[indexF], indexE);

        VertexType *testV=mesh.face[indexF].V(indexE);
        std::vector<FaceType*> face_vec;
        std::vector<int> edge_vec;

        vcg::face::VFOrderedStarFF(start_pos,face_vec,edge_vec);

        for (size_t i=0;i<face_vec.size();i++)
        {
            assert((face_vec[i]->V(0)==testV)||
                   (face_vec[i]->V(1)==testV)||
                   (face_vec[i]->V(2)==testV));

            //get the edge that doesnt' have such vertex
            int indexE=1;
            if (face_vec[i]->V(1)==testV)indexE=2;
            if (face_vec[i]->V(2)==testV)indexE=0;
            vcg::face::Pos<FaceType> EdgeP(face_vec[i],indexE);
            //check if it is border
            if (EdgeP.IsBorder())continue;
            EdgeP.FlipF();
            surroundEdges.push_back(EdgeP);
        }
    }

    /**
     * @brief getNodesSurroundingSing returns the outside nodes around a given singular vertex node.
     * @param indexF
     * @param indexE
     * @param surroundNodes
     */
    void getNodesSurroundingSing(const vcg::face::Pos<FaceType> &singPos,
                                 std::vector<std::vector<size_t> > &surroundNodes,
                                 std::vector<std::vector<size_t> > &surroundDir,
                                 std::vector<CoordType> &surroundPos,
                                 std::vector<size_t> &surroundFace,
                                 std::vector<size_t> &surroundEdge)
    {
        surroundNodes.clear();
        surroundDir.clear();
        surroundPos.clear();
        surroundFace.clear();
        surroundEdge.clear();

        // get the surrounding edges
        std::list<vcg::face::Pos<FaceType> > surrPos;
        getPosSurroundingESing(singPos, surrPos);

        //then save for each coordinate the nodes

        // get the corresponding nodes
        std::map<CoordType, NodeList> surroundCoordNodes;
        faceEdgeNodes(surrPos, surroundCoordNodes);

        //transform into a vector
        typename std::map<CoordType, NodeList>::iterator IteM;
        for (IteM=surroundCoordNodes.begin();IteM!=surroundCoordNodes.end();IteM++)
        {
            surroundNodes.push_back(std::vector<size_t>());
            surroundDir.push_back(std::vector<size_t>());

            bool addedFirst=false;
            for (NodeListIterator NIte=(*IteM).second.begin();NIte!=(*IteM).second.end();NIte++)
            {
                size_t IndexN=(*NIte);
                surroundNodes.back().push_back(IndexN);

                assert(IsEdgeNode(IndexN));
                EdgeGeodesicNode *edgeN = static_cast<EdgeGeodesicNode*>(Graph[IndexN]);

                //get face and edge index
                size_t IndexF=edgeN->faceID;
                size_t IndexE=edgeN->edgeIdx;
                size_t M4Dir=edgeN->M4Dir;
                CoordType Pos=edgeN->pos;
                surroundDir.back().push_back(M4Dir);

                if (!addedFirst)
                {
                    surroundPos.push_back(Pos);
                    surroundFace.push_back(IndexF);
                    surroundEdge.push_back(IndexE);
                    addedFirst=true;
                }
                else //check should be the same face and edge and pos for all nodes
                {
                    assert(Pos==surroundPos.back());
                    assert(IndexF==surroundFace.back());
                    assert(IndexE==surroundEdge.back());
                }
                //surroundNodes.push_back(std::vector<size_t>((*IteM).second.begin(),(*IteM).second.end()));
            }
        }
    }



    void InterpolateNodesUV(const std::vector<CoordType> &NodePos,
                            const std::vector<size_t> &surroundFace,
                            const std::vector<size_t> &surroundEdge,
                            std::vector<vcg::Point2<ScalarType> > &surroundNodesUV)
    {
        assert(NodePos.size()==surroundFace.size());
        assert(surroundFace.size()==surroundEdge.size());

        surroundNodesUV.clear();
        for (size_t i=0;i<NodePos.size();i++)
        {
            //get face and edge index
            size_t IndexF=surroundFace[i];
            size_t IndexE=surroundEdge[i];

            //then interpolate the UV coords
            VertexType *v0=mesh.face[IndexF].V0(IndexE);
            VertexType *v1=mesh.face[IndexF].V1(IndexE);

            //get the interpolation value
            ScalarType L_edge=(v0->P()-v1->P()).Norm();
            ScalarType L_0=(v0->P()-NodePos[i]).Norm();
            ScalarType alpha=1-L_0/L_edge;
            assert(alpha>=0);
            assert(alpha<=1);

            //interpolate the UV coordinates
            vcg::Point2<ScalarType> InterpUV=v0->T().P()*alpha+v1->T().P()*(1-alpha);

            //normalize the direction
            InterpUV.Normalize();

            //get the two UV coords
            surroundNodesUV.push_back(InterpUV);
        }
    }

    void GetSingDirections(const std::vector<size_t> &SingIndex,
                           std::vector<vcg::Point2<ScalarType> > &UVSingDirection,
                           std::vector<size_t> &SingExitFace,
                           std::vector<CoordType> &SingDirections)
    {

        //then convert each sing to UV
        UVSingDirection.clear();
        SingExitFace.clear();
        SingDirections.clear();

        //then cycle over all the singularities
        for (size_t i=0;i<SingIndex.size();i++)
        {
            int SingNode=SingIndex[i];
            assert(IsSingularity(SingNode));
            SingularityGeodesicNode *singN = static_cast<SingularityGeodesicNode*>(Graph[SingNode]);

            //set the origin
            CoordType Origin=singN->pos;

            CoordType SingDir=singN->mainDirection;
            SingDir.Normalize();
            size_t exitFace=singN->exitFace;

            //create a triangle of the face
            FaceType *ExitF=&mesh.face[exitFace];

            //push per sign indexes
            SingDirections.push_back(SingDir);
            SingExitFace.push_back(exitFace);

            //then create the triangle of the face
            vcg::Triangle3<ScalarType> InterpFace(ExitF->P(0)-Origin,
                                                  ExitF->P(1)-Origin,
                                                  ExitF->P(2)-Origin);

            //get the triangle in UV Space
            vcg::Triangle2<ScalarType> UVFace(ExitF->V(0)->T().P(),
                                              ExitF->V(1)->T().P(),
                                              ExitF->V(2)->T().P());

            //scale the direction by the barycenter
            CoordType Bary=(InterpFace.P(0)+InterpFace.P(1)+InterpFace.P(2))/3;
            SingDir*=Bary.Norm();

            //find the barycentric coords
            CoordType L;
            vcg::InterpolationParameters(InterpFace,SingDir,L);

            //then interpolate UV
            vcg::Point2<ScalarType> UVDir=L.X()*UVFace.P(0)+
                    L.Y()*UVFace.P(1)+
                    L.Z()*UVFace.P(2);
            UVDir.Normalize();
            UVSingDirection.push_back(UVDir);
        }
    }

    void AssociateUVPosUVDirections(const std::vector<vcg::Point2<ScalarType> > &surroundNodesUV,
                                    const std::vector<vcg::Point2<ScalarType> > &UVSingDirection,
                                    std::vector<int> &BestD)
    {
        //initialization of the vector
        BestD=std::vector<int>(surroundNodesUV.size(),-1);
        //then cycle over all surronding nodes
        for (size_t i=0;i<surroundNodesUV.size();i++)
        {
            //find the most coherent direction
            ScalarType BestProd=-1;
            vcg::Point2<ScalarType> Dir0=surroundNodesUV[i];
            //check versus all sing
            assert(UVSingDirection.size()>0);
            for (size_t j=0;j<UVSingDirection.size();j++)
            {
                vcg::Point2<ScalarType> Dir1=UVSingDirection[j];
                ScalarType TestProd=(Dir0*Dir1);
                if (TestProd<BestProd)continue;
                BestProd=TestProd;
                BestD[i]=j;
            }
            assert(BestD[i]!=-1);
        }
    }

    void FindExpectedDir(const std::vector<int> &BestD,
                         const std::vector<CoordType> &SingDirections,
                         const std::vector<size_t> &SingExitFace,
                         const std::vector<size_t> &surroundFace,
                         std::vector<size_t> &ExpectedDir)
    {
        for (size_t i=0;i<BestD.size();i++)
        {
            size_t indexSingN=BestD[i];
            assert(indexSingN>=0);
            assert(indexSingN<SingDirections.size());
            assert(indexSingN<SingExitFace.size());

            //get the face of the node
            size_t IndexF=surroundFace[i];
            assert(IndexF>=0);
            assert(IndexF<mesh.face.size());

            FaceType *FaceNode=&mesh.face[IndexF];

            //then get the expected direction
            assert(indexSingN<SingDirections.size());
            assert(indexSingN<SingExitFace.size());
            assert(indexSingN>=0);

            //then get the directions of the singularity
            CoordType SingD=SingDirections[indexSingN];
            IndexF=SingExitFace[indexSingN];
            assert(IndexF>=0);
            assert(IndexF<mesh.face.size());
            FaceType *FaceSing=&mesh.face[IndexF];

            //get the direction in the target face
            size_t ExpDir=vcg::tri::CrossField<MeshType>::FollowDirectionI(*FaceSing,*FaceNode,SingD);

            //push the expected direction
            ExpectedDir.push_back(ExpDir);
        }
    }

    void AttachSingToNodes(const std::vector<size_t> &SingIndex,
                           const std::vector<std::vector<size_t> > &surroundNodes,
                           const std::vector<int> &BestD,
                           const std::vector<size_t> &ExpectedDir,
                           const std::vector<std::vector<size_t> > &surroundDir)
    {
        assert(surroundNodes.size()==BestD.size());
        assert(surroundNodes.size()==ExpectedDir.size());
        assert(surroundNodes.size()==surroundDir.size());

        for (size_t i=0;i<surroundNodes.size();i++)
        {
            assert(BestD[i]>=0);
            assert(BestD[i]<(int)SingIndex.size());
            assert(surroundNodes[i].size()==surroundDir[i].size());

            size_t SingN=SingIndex[BestD[i]];

            //then find the node which have the same direction
            for (size_t j=0;j<surroundNodes[i].size();j++)
            {
                if (surroundDir[i][j]!=ExpectedDir[i])continue;

                int SurrN=surroundNodes[i][j];

                CoordType posSing=NodePos(SingN);
                CoordType posNode=NodePos(SurrN);

                CoordType TestDir=posNode-posSing;
                ScalarType length = TestDir.Norm();
                TestDir.Normalize();

                CoordType Direction=getDirection(SingN);
                ScalarType angle = vcg::Angle(Direction, TestDir) * 180 / M_PI;
                //            //if (angle1 <= ScalarType(90) && angle2 <= /*max_angle*/ScalarType(90))
                Graph[SingN]->neighbors.push_back(NeighborInfo(SurrN, angle, length));
            }
        }
    }

    void ConnectSing(const vcg::face::Pos<FaceType> &SingP,
                     const std::vector<size_t> &SingIndex)
    {
        //parametrize the star of faces
        ParametrizeStar(SingP);

        //safety check
        assert(SingP.V()->T().P()==vcg::Point2<ScalarType>(0,0));

        //get the surrounding nodes
        std::vector<std::vector<size_t> > surroundNodes;
        std::vector<std::vector<size_t> > surroundDir;
        std::vector<CoordType> surroundPos;
        std::vector<size_t> surroundFace;
        std::vector<size_t> surroundEdge;

        getNodesSurroundingSing(SingP,surroundNodes,surroundDir,surroundPos,surroundFace,surroundEdge);
        assert(surroundNodes.size()==surroundDir.size());
        assert(surroundNodes.size()==surroundPos.size());
        assert(surroundNodes.size()==surroundFace.size());
        assert(surroundNodes.size()==surroundEdge.size());

        //then find the UV coords per Node considering they have the
        //same position only one per coord is needed
        std::vector<vcg::Point2<ScalarType> > surroundNodesUV;
        InterpolateNodesUV(surroundPos,surroundFace,surroundEdge,surroundNodesUV);
        assert(surroundNodesUV.size()==surroundNodes.size());

        //then find singularities directions in UV
        std::vector<vcg::Point2<ScalarType> > UVSingDirection;
        std::vector<size_t> SingExitFace;
        std::vector<CoordType> SingDirections;
        GetSingDirections(SingIndex,UVSingDirection,SingExitFace,SingDirections);
        assert(UVSingDirection.size()==SingIndex.size());
        assert(SingExitFace.size()==SingIndex.size());
        assert(SingDirections.size()==SingIndex.size());

        //associate for each node coord the closest direction in UV space
        std::vector<int> BestD;
        AssociateUVPosUVDirections(surroundNodesUV,UVSingDirection,BestD);
        assert(BestD.size()==surroundNodesUV.size());

        //finally find the direction which is congruent with the separatrix
        std::vector<size_t> ExpectedDir;

        //then find the expected following direction
        FindExpectedDir(BestD,SingDirections,SingExitFace,surroundFace,ExpectedDir);

        //finally add the links
        AttachSingToNodes(SingIndex,surroundNodes,BestD,ExpectedDir,surroundDir);
    }

    /**
     * @brief addSingularities adds singularity nodes and sets their connections.
     */
    void addSingularities() {
        // query if an attribute is present or not
        assert(vcg::tri::HasPerVertexAttribute(mesh,std::string("Singular")));
        assert(vcg::tri::HasPerVertexAttribute(mesh,std::string("SingularIndex")));
        typename MeshType::template PerVertexAttributeHandle<bool> Handle_Singular;
        Handle_Singular = vcg::tri::Allocator<MeshType>::template GetPerVertexAttribute<bool>(mesh,std::string("Singular"));
        typename MeshType::template PerVertexAttributeHandle<int> Handle_SingularIndex;
        Handle_SingularIndex = vcg::tri::Allocator<MeshType>::template GetPerVertexAttribute<int>(mesh,std::string("SingularIndex"));

        PerVertNodes.clear();
        PerVertNodes.resize(mesh.vert.size());
        vcg::tri::UpdateFlags<MeshType>::VertexClearS(mesh);

        vcg::face::Pos<FaceType> currP;
        for (size_t i = 0; i < mesh.face.size(); ++i)
            for (int j = 0; j < mesh.face[i].VN(); ++j) {
                // already explored
                if (mesh.face[i].V(j)->IsS())
                    continue;
                // select to set as explored
                mesh.face[i].V(j)->SetS();
                // not a singular
                if (!Handle_Singular[*mesh.face[i].V(j)])
                    continue;

                //check non border
                assert(!mesh.face[i].V(j)->IsB());

                // find the separatrices
                currP.Set(&mesh.face[i], j, mesh.face[i].V(j));
                //int valence = mesh.expectedValence(*mesh.face[i].V(j));
                int valence= vcg::tri::CrossField<MeshType>::expectedValence(mesh,*mesh.face[i].cV(j));
                std::vector<FaceType*> Faces;
                std::vector<CoordType> directions;

                size_t sampledDir=vcg::tri::CrossField<MeshType>::FindSeparatrices(currP,directions,Faces,WrongTris,valence);

                if (sampledDir!=directions.size())
                    std::cout << "Sampled " << sampledDir << "directions instead of " << valence << " CORRECTED!" <<std::endl;

                // add the nodes
                std::vector<size_t> newSingNodes;
                addSingularityNodes(vcg::tri::Index(mesh, mesh.face[i].V(j)), directions,Faces, newSingNodes);


                //if it is not Ok then insert into the ones needed to be visualized
                if (sampledDir!=directions.size()){
                    WrongNodes.insert(WrongNodes.end(),newSingNodes.begin(),newSingNodes.end());
                }

                //conect singularities with the rest of the nodes
                ConnectSing(currP,newSingNodes);
            }
        vcg::tri::UpdateFlags<MeshType>::VertexClearS(mesh);
    }

    /**
     * @brief resetPerEdgeNodes collects indices of nodes for each edge
     */
    void resetPerEdgeNodes() {
        // allocate space
        PerEdgeNodes.clear();
        PerEdgeNodes.resize(mesh.face.size());
        for (size_t i = 0; i < PerEdgeNodes.size(); ++i)
            PerEdgeNodes[i].resize(mesh.face[i].VN(), std::vector<NodeList>(4));
        // store the indices
        for (size_t i = 0; i < nEdgeNodes; ++i) {
            assert(Graph[i]->nodeType == EdgeNodeType);
            EdgeGeodesicNode *edgeN = static_cast<EdgeGeodesicNode*>(Graph[i]);
            PerEdgeNodes[edgeN->faceID][edgeN->edgeIdx][edgeN->M4Dir].push_back(i);
        }
    }


public:

    void GetUnreferencedNodes(std::vector<int> &Unreferenced)
    {
        std::vector<size_t> referenceCount(NumNodes(),0);
        for (size_t i = 0; i < Graph.size(); ++i)
            for (NeighborsIterator nIt = Graph[i]->neighbors.begin(); nIt != Graph[i]->neighbors.end(); ++nIt)
            {
                if (!nIt->Valid)continue;

                assert(nIt->nodeID < referenceCount.size());
                assert(nIt->nodeID >=0);

                ++referenceCount[nIt->nodeID];
            }

        for (size_t i = 0; i < Graph.size(); ++i)
        {
            if (referenceCount[i]>0)continue;
            Unreferenced.push_back(i);
        }
    }

    bool IsDeadEnd(int IndexNode)
    {
        assert(IndexNode>=0);
        assert(IndexNode<Graph.size());
        size_t NumArcs=0;
        for (NeighborsIterator nIt = Graph[IndexNode]->neighbors.begin();
             nIt != Graph[IndexNode]->neighbors.end(); ++nIt)
        {
            if (!nIt->Valid)continue;
            NumArcs++;
        }
        if (NumArcs>0)return false;
        return true;
    }

    void GetDeadEndNodes(std::vector<int> &DeadEnd)
    {
        for (size_t i = 0; i < Graph.size(); ++i)
        {
            if (IsDeadEnd(i))
                DeadEnd.push_back(i);
        }
    }

private:
    /**
     * @brief deleteUnreferencedNodes discards edge nodes that are not reached by any other node.
     * @param erase_dead_end
     */
    void deleteUnreferencedNodes(const bool erase_dead_end = true) {


        std::vector<size_t> referenceCount, newIndexes;
        std::list<BaseGeodesicNode*> swap;
        size_t newNEdgeNodes = 0;
        bool modified=false;
        do {
            std::cout<<"Deleting Unreferenced Step"<<std::endl;
            modified=false;
            // count the references
            referenceCount.assign(Graph.size(), 0);
            for (size_t i = 0; i < Graph.size(); ++i)
                for (NeighborsIterator nIt = Graph[i]->neighbors.begin(); nIt != Graph[i]->neighbors.end(); ++nIt)
                    if (nIt->nodeID < referenceCount.size())
                        ++referenceCount[nIt->nodeID];

            // swap new nodes and save the new index of vertices
            newIndexes.assign(Graph.size(), std::numeric_limits<size_t>::max());
            swap.clear();
            newNEdgeNodes = 0;
            for (size_t i = 0; i < nEdgeNodes; ++i)
            {
                if (!Graph[i])
                    continue;
                assert(Graph[i]->nodeType == EdgeNodeType);
                EdgeGeodesicNode *node=static_cast<EdgeGeodesicNode*>(Graph[i]);
                bool IsBorderNode=node->IsBorder;
                //delete the ones that have no reference or dead end
                //the ones on the border are kept!
                if ((!IsBorderNode)&&
                        ((referenceCount[i] == 0) ||
                         (erase_dead_end && Graph[i]->neighbors.empty())))
                {
                    delete Graph[i];
                    Graph[i] = 0;
                    modified=true;
                    continue;
                }
                swap.push_back(Graph[i]);
                newIndexes[i] = swap.size() - 1;
                ++newNEdgeNodes;
            }
            //save the index of remaining vertices
            for (size_t i = nEdgeNodes; i < Graph.size(); ++i)
            {
                assert(Graph[i]->nodeType != EdgeNodeType);
                swap.push_back(Graph[i]);
                newIndexes[i] = swap.size() - 1;
            }

            // delete references to deleted nodes
            //if (erase_dead_end)
            for (size_t i = 0; i < Graph.size(); ++i)
                if (Graph[i])//if it has not been deleted
                    for (NeighborsIterator nIt = Graph[i]->neighbors.begin(); nIt != Graph[i]->neighbors.end(); )
                    {
                        //delete null neighbors
                        if (!Graph[nIt->nodeID]) /*||
                                    (Graph[nIt->nodeID]->nodeType == EdgeNodeType && Graph[nIt->nodeID]->neighbors.empty()))*/
                            nIt = Graph[i]->neighbors.erase(nIt);
                        else
                            ++nIt;
                    }

            // swap buffers and updated references
            Graph.assign(swap.begin(), swap.end());
            nEdgeNodes = newNEdgeNodes;
            for (size_t i = 0; i < Graph.size(); ++i)
            {
                //this should has been already deleted
                //assert(newIndexes[i]!=std::numeric_limits<size_t>::max());
                for (NeighborsIterator nIt = Graph[i]->neighbors.begin(); nIt != Graph[i]->neighbors.end(); )
                    //                    if (newIndexes[nIt->nodeID] >= Graph.size())
                    //                        nIt = Graph[i]->neighbors.erase(nIt);
                    //                    else {
                {
                    //this should has been already deleted
                    assert(newIndexes[nIt->nodeID]!=std::numeric_limits<size_t>::max());
                    nIt->nodeID = newIndexes[nIt->nodeID];
                    ++nIt;
                }
            }

            //update opposite node of sigularities ans sinks
            for (size_t i = 0; i < Graph.size(); ++i)
            {
                //chenge the opposite of the same face
                if (IsEdgeNode(i))//continue;
                {
                    //get the sink node
                    EdgeGeodesicNode *currEdgeNode=static_cast<EdgeGeodesicNode*>(Graph[i]);

                    //get the relative singularity
                    int oppNodeSameFace=currEdgeNode->OppositeSameFID;
                    if (oppNodeSameFace>=0)
                    {
                        if (newIndexes[oppNodeSameFace]==std::numeric_limits<size_t>::max())
                            currEdgeNode->OppositeSameFID=-1;
                        else
                            currEdgeNode->OppositeSameFID=newIndexes[oppNodeSameFace];
                    }
                }
                else
                    if (IsSink(i))
                    {
                        //get the sink node
                        SinkGeodesicNode *currSink=static_cast<SinkGeodesicNode*>(Graph[i]);
                        //get the relative singularity
                        size_t oppNode=currSink->oppositeSingularityNode;
                        currSink->oppositeSingularityNode=newIndexes[oppNode];
                    }else
                    {
                        assert(IsSingularity(i));
                        //get the sink node
                        SingularityGeodesicNode *currSing=static_cast<SingularityGeodesicNode*>(Graph[i]);

                        //get the relative singularity
                        size_t oppNode=currSing->oppositeSinkNode;
                        currSing->oppositeSinkNode=newIndexes[oppNode];
                    }
            }
            //update wrong nodes
            for (size_t i=0;i<WrongNodes.size();i++)
                WrongNodes[i]=newIndexes[WrongNodes[i]];
            for (size_t i=0;i<NoOppNodes.size();i++)
                NoOppNodes[i]=newIndexes[NoOppNodes[i]];
            // update the PerEdgeNodes table
            //resetPerEdgeNodes();
        } while (modified);

        resetPerEdgeNodes();
    }

public:

    /**
     * @brief getDirection returns the direction of a given node.
     * @param indexNode
     * @return
     */
    CoordType getDirection(size_t indexNode) {
        assert(indexNode < Graph.size());

        EdgeGeodesicNode *edgeN = 0;
        SingularityGeodesicNode *singN = 0;
        SinkGeodesicNode *sinkN = 0;
        switch (Graph[indexNode]->nodeType) {
        case EdgeNodeType:
            edgeN = static_cast<EdgeGeodesicNode*>(Graph[indexNode]);
            return vcg::tri::CrossField<MeshType>::CrossVector(mesh.face[edgeN->faceID], edgeN->M4Dir);

        case SinkNodeType:
            sinkN = static_cast<SinkGeodesicNode*>(Graph[indexNode]);
            indexNode = sinkN->oppositeSingularityNode;
            singN = static_cast<SingularityGeodesicNode*>(Graph[indexNode]);
            return (-singN->mainDirection);

        case SingularityNodeType:
            singN = static_cast<SingularityGeodesicNode*>(Graph[indexNode]);
            return singN->mainDirection;

        default:
            return CoordType(0, 0 ,0);
        }
    }

    void SelectUnreferencedNodes()
    {
        for (size_t i=0;i<Graph.size();i++)
            Graph[i]->Selected=false;

        std::vector<size_t> referenceCount;

        // count the references
        referenceCount.assign(Graph.size(), 0);
        for (size_t i = 0; i < Graph.size(); ++i)
            for (NeighborsIterator nIt = Graph[i]->neighbors.begin(); nIt != Graph[i]->neighbors.end(); ++nIt)
                ++referenceCount[nIt->nodeID];

        for (size_t i = 0; i < nEdgeNodes; ++i)
        {
            assert(Graph[i]->nodeType == EdgeNodeType);
            Graph[i]->Selected=(referenceCount[i] == 0 );
        }
    }

    void SelectDeadEndNodes()
    {
        for (size_t i=0;i<Graph.size();i++)
            Graph[i]->Selected=false;

        for (size_t i = 0; i < nEdgeNodes; ++i)
        {
            assert(Graph[i]->nodeType == EdgeNodeType);
            Graph[i]->Selected=Graph[i]->neighbors.empty();
        }
    }
private:
    /**
     * @brief addSinks adds a sink node for each singularity node.
     */
    void addSinks() {
        for (size_t i = nEdgeNodes; i < nEdgeNodes + nSingNodes; ++i) {
            assert(Graph[i]->nodeType == SingularityNodeType);
            SingularityGeodesicNode *SingNode = static_cast<SingularityGeodesicNode*>(Graph[i]);

            // new node
            Graph.push_back(new SinkGeodesicNode(i, SingNode->pos));

            SingNode->oppositeSinkNode=Graph.size()-1;

            //Graph.back()->oppositeSingularityNode=SingNode->
            PerVertNodes[SingNode->vertexID].push_back(Graph.size() - 1);
            // set the connections
            for (NeighborsIterator nIt = SingNode->neighbors.begin(); nIt != SingNode->neighbors.end(); ++nIt) {
                assert(Graph[nIt->nodeID]->nodeType == EdgeNodeType);
                int oppNode = ComputeOppositeEdgeNode(nIt->nodeID,2);//oppositeNode(nIt->nodeID, -getDirection(nIt->nodeID));
                assert(oppNode != -1);
                // link to the sink
                Graph[oppNode]->neighbors.push_back(NeighborInfo(Graph.size() - 1, nIt->angle, nIt->length));
            }
        }
    }

    void UpdateStats()
    {
        //cumulate all nodes info
        Stats.NumNodes=Graph.size();
        Stats.CorrectedSeparatrix=WrongNodes.size();
        Stats.WrongOpposite=NoOppNodes.size();
        Stats.NumSingularities=0;
        for (size_t i=0;i<mesh.vert.size();i++)
        {
            //if (!mesh.IsSingularVert(mesh.vert[i]))continue;
            if (vcg::tri::CrossField<MeshType>::IsSingular(mesh,mesh.vert[i]))continue;
            Stats.NumSingularities++;
        }

        Stats.EdgeNodes=0;
        Stats.BoundaryNodes=0;
        Stats.SingularityNodes=0;
        Stats.SinkNodes=0;

        for (size_t i=0;i<Graph.size();i++)
        {
            if (IsEdgeNode(i))Stats.EdgeNodes++;
            if (IsSingularity(i))Stats.SingularityNodes++;
            if (IsSink(i))Stats.SinkNodes++;
            if (IsEdgeNode(i)&&IsBorderNode(i))Stats.BoundaryNodes++;
        }
    }

public:
    /**
     * @brief InitGraph builds the graph by sampling the edge nodes, adding singularity and sink nodes and connections among all.
     * @param subDivision
     * @param maxDev
     */
    void InitGraph(const int subDivision = 5,
                   const ScalarType maxDev = ScalarType(45),
                   bool AddSingNodes=true,
                   bool RemoveUnref=true)
    {



        //update border flag
        int t_init=clock();

        WrongNodes.clear();
        NoOppNodes.clear();

        if (AddSingNodes)
            SplitAdjSing();

        vcg::tri::UpdateFlags<MeshType>::FaceBorderFromFF(mesh);
        vcg::tri::UpdateFlags<MeshType>::VertexBorderFromFaceBorder(mesh);

        //SplitSing();

        //SplitSing();
        //SplitSing();
        //while (SplitSing());

        max_angle = maxDev;

        // compute the target sample distance
        ScalarType minSize, maxSize, avgSize;
        edgeSize(minSize, maxSize, avgSize);
        ScalarType sampleDistance = avgSize / static_cast<ScalarType>(subDivision);

        // initialize the graph
        clearGraph();

        // sample edge nodes
        int t0 = clock();
        std::cout << "Initializing graph ... Subdivision "<< subDivision <<
                     " Edge distance "<<sampleDistance<< std::endl;
        std::cout.flush();
        sampleMesh(subDivision,sampleDistance);
        int t1 = clock();
        std::cout << "done. Sampled " << nEdgeNodes << " edge-nodes. Time " << t1 - t0 << std::endl;

        // set connections
        t0 = clock();
        std::cout << "Setting connections ... ";
        std::cout.flush();
        initConnections();
        t1 = clock();
        std::cout << "done. Time " << t1 - t0 << std::endl;

        std::cout << "Initializing points-nodes mapping ... ";
        std::cout.flush();
        t0 = clock();
        initCoordTable();
        t1 = clock();
        std::cout << "done. Time " << t1 - t0 << std::endl;

        std::cout << "Adding singularity nodes ... ";
        std::cout.flush();
        t0 = clock();
        if (AddSingNodes)
            addSingularities();
        t1 = clock();
        std::cout << "done. Added " << nSingNodes << " singularity nodes. Time " << t1 - t0 << std::endl;

        std::cout << "Adding sink nodes ...";
        std::cout.flush();
        t0 = clock();
        addSinks();
        t1 = clock();
        std::cout << "done. Time " << t1 - t0 << std::endl;

        if (RemoveUnref)
        {
            std::cout << "Deleting unconnected nodes ...";
            std::cout.flush();
            t0 = clock();
            size_t size0 = Graph.size();
            deleteUnreferencedNodes();
            t1 = clock();
            std::cout << "done. Deleted " << size0 - Graph.size() << " nodes. Time " << t1 - t0 << std::endl;
        }

        std::cout << "Re-initializing points-nodes mapping ... ";
        std::cout.flush();
        t0 = clock();
        initCoordTable();
        t1 = clock();
        std::cout << "done. Time " << t1 - t0 << std::endl;

        std::cout << "Initialize  opposite nodes... ";
        std::cout.flush();
        t0 = clock();
        InitOppositeNodes();
        t1 = clock();
        std::cout << "done. Time " << t1 - t0 << std::endl;

        int t_final=clock();
        //update stats
        UpdateStats();
        Stats.timeGraph=(ScalarType)(t_final-t_init)/(ScalarType)(CLOCKS_PER_SEC);
    }

    //private:
    std::vector<size_t>         Father;      ///< Vector of indices: for each node, the index of its father in geodesic distance computed from the current seed node.
    std::vector<size_t>         Source;      ///< Vector of indices: for each node, the index of the current seed node.
    std::vector<ScalarType>     Distance;    ///< Vector of distances of each node from the current seed.
    std::vector<size_t>         Group;       ///< Vector of indices: for each node the group of sources it belongs.



public:



    /**
     * @brief GetSingularitiesInfo returns all info on singularities
     * @param SingNode
     * @param SinkNode
     * @param meshVert
     * @param expValence
     */
    void GetSingularitiesInfo(std::vector<std::vector<size_t > > &SingNode,
                              std::vector<std::vector<size_t > > &SinkNode,
                              std::vector<size_t> &meshVert,
                              std::vector<size_t> &expValence)
    {
        SingNode.clear();
        SinkNode.clear();
        meshVert.clear();
        expValence.clear();

        std::map<int,int> SingMap;

        //then run over all the sinks nodes
        for (size_t i = nEdgeNodes + nSingNodes; i < Graph.size(); ++i)
        {
            assert(Graph[i]->nodeType == SinkNodeType);

            //get the sink node
            SinkGeodesicNode *currSink=(SinkGeodesicNode *)Graph[i];

            //get the relative singularity
            size_t oppositeSingularityNode=currSink->oppositeSingularityNode;
            assert((oppositeSingularityNode>=0)&&(oppositeSingularityNode<Graph.size()));
            assert(Graph[oppositeSingularityNode]->nodeType == SingularityNodeType);

            //get the singularity Node
            SingularityGeodesicNode *currSing=(SingularityGeodesicNode *)Graph[oppositeSingularityNode];

            //get the index of the vertex
            size_t SingVertexID=currSing->vertexID;
            assert(SingVertexID>=0);
            assert(SingVertexID<mesh.vert.size());

            //and get the expected valence
            //int expectedSingValue=mesh.expectedValence(mesh.vert[SingVertexID]);//expectedValence(mesh.vert[SingVertexID]);
            int expectedSingValue=vcg::tri::CrossField<MeshType>::expectedValence(mesh,mesh.vert[SingVertexID]);//expectedValence(mesh.vert[SingVertexID]);

            //then get the index if it exist
            int curr_index=-1;
            if (SingMap.count(SingVertexID)==0)
            {
                //add the entries
                meshVert.push_back(SingVertexID);
                expValence.push_back(expectedSingValue);
                curr_index=meshVert.size()-1;

                //then create the map
                SingMap[SingVertexID]=curr_index;

                //add a new entry
                SingNode.push_back(std::vector<size_t >());
                SinkNode.push_back(std::vector<size_t >());
                assert((int)SingNode.size()==curr_index+1);
                assert((int)SinkNode.size()==curr_index+1);

            }
            else
                curr_index=SingMap[SingVertexID];

            assert(SinkNode.size()==SingNode.size());
            assert((curr_index>=0)&&
                   (curr_index<(int)SinkNode.size())&&
                   (curr_index<(int)meshVert.size())&&
                   (curr_index<(int)expValence.size()));

            SingNode[curr_index].push_back(oppositeSingularityNode);
            SinkNode[curr_index].push_back(i);
        }
    }

    void GetSingularityNodes(std::vector<size_t > &SingNode,
                             std::vector<size_t > &SinkNode)
    {
        SingNode.clear();
        SinkNode.clear();
        std::vector<std::vector<size_t > > SingNodeSwap,SinkNodeSwap;
        std::vector<size_t > meshVert,expValence;
        GetSingularitiesInfo(SingNodeSwap,SinkNodeSwap,meshVert,expValence);

        //then put in the vector
        for (size_t i=0;i<SingNodeSwap.size();i++)
        {
            SingNode.insert(SingNode.end(),SingNodeSwap[i].begin(),SingNodeSwap[i].end());
            SinkNode.insert(SinkNode.end(),SinkNodeSwap[i].begin(),SinkNodeSwap[i].end());
        }
    }


#ifndef NO_TRACING_OPENGL

    // drawing functions
    void GlDrawWrongTris()
    {
        glPushAttrib(GL_ALL_ATTRIB_BITS);
        glDisable(GL_LIGHTING);
        glDepthRange(0,0.99999);
        glLineWidth(3);
        for (size_t i=0;i<WrongTris.size();++i)
        {
            glBegin(GL_LINE_LOOP);
            vcg::glColor(vcg::Color4b(255,0,255,255));
            vcg::glVertex(WrongTris[i].P(0));
            vcg::glVertex(WrongTris[i].P(1));
            vcg::glVertex(WrongTris[i].P(2));
            glEnd();
        }

        glEnable(GL_LIGHTING);
        for (size_t i=0;i<WrongNodes.size();++i)
            GlDrawNodeDebugInfo(WrongNodes[i],45,100,0.001,true,false);

        for (size_t i=0;i<NoOppNodes.size();++i)
            GlDrawNodeDebugInfo(NoOppNodes[i],45,100,0.001,true,false,vcg::Color4b(255,0,0,255));

        glPopAttrib();
    }

    //    // drawing functions
    //    void GlDrawPath(const Path &Path,
    //                    vcg::Color4b path_col=vcg::Color4b(255,0,0,255),
    //                    ScalarType Line_Width=10,
    //                    bool loop=false)
    //    {
    //        vcg::glColor(path_col);
    //        glPushAttrib(GL_ALL_ATTRIB_BITS);
    //        glDisable(GL_LIGHTING);
    //        glDepthRange(0,0.99999);
    //        glLineWidth(Line_Width);
    //        glBegin(GL_LINE_STRIP);
    //        for (size_t i=0;i<Path.nodes.size();++i)
    //        {
    //            //        if (i==0)
    //            //            vcg::glColor(vcg::Color4b(255,0,0,255));
    //            //        else
    //            //
    //            vcg::glColor(path_col);
    //            vcg::glVertex(Graph[Path.nodes[i]]->pos);
    //        }
    //        if (loop)
    //            vcg::glVertex(Graph[Path.nodes[0]]->pos);
    //        glEnd();
    //        glPopAttrib();
    //    }


    // drawing functions
    void GlDrawPath(const Path &Path,
                    vcg::Color4b path_col=vcg::Color4b(255,0,0,255),
                    ScalarType Line_Width=10,
                    bool loop=false)
    //bool black_bound=false)
    {
        vcg::glColor(path_col);
        glPushAttrib(GL_ALL_ATTRIB_BITS);
        glDisable(GL_LIGHTING);
        glDepthRange(0,0.99999);
        glLineWidth(Line_Width);
        glBegin(GL_LINE_STRIP);
        for (size_t i=0;i<Path.nodes.size();++i)
        {
            vcg::glColor(path_col);
            vcg::glVertex(Graph[Path.nodes[i]]->pos);
        }
        if (loop)
            vcg::glVertex(Graph[Path.nodes[0]]->pos);
        glEnd();
        //        if (black_bound)
        //        {
        //            glDepthRange(0,0.999999);
        //            glLineWidth(Line_Width+4);
        //            glBegin(GL_LINE_STRIP);
        //            for (size_t i=0;i<Path.nodes.size();++i)
        //            {
        //                vcg::glColor(vcg::Color4b(0,0,0,255));
        //                vcg::glVertex(Graph[Path.nodes[i]]->pos);
        //            }
        //            if (loop)
        //                vcg::glVertex(Graph[Path.nodes[0]]->pos);
        //            glEnd();
        //        }
        glPopAttrib();
    }

    void GlDrawPaths(const std::vector<Path> &Paths)
    {
        for (size_t i=0;i<Paths.size();++i)
        {
            vcg::Color4b path_col=vcg::Color4b::Scatter(Paths.size(),i);
            GlDrawPath(Paths[i],path_col);
        }
    }

    void GlDrawConnections(size_t Node,
                           const ScalarType &max_deviation,
                           const ScalarType &penalty_drift)
    {
        glPushAttrib(GL_ALL_ATTRIB_BITS);
        glDisable(GL_LIGHTING);
        glDepthRange(0,0.99999);

        glLineWidth(2);
        vcg::Point3f pos0;
        pos0.Import<ScalarType>(Graph[Node]->pos);

        ScalarType MaxGeoL,MinGeoL;

        Graph[Node]->MinMaxDistance(MinGeoL,MaxGeoL,max_deviation,penalty_drift);
        //Graph[Node]->MinMaxAngle(MinGeoL,MaxGeoL);

        glBegin(GL_LINES);
        for (NeighborsIterator nIt = Graph[Node]->neighbors.begin(); nIt != Graph[Node]->neighbors.end(); ++nIt)
        {
            if (!nIt->Valid)continue;


            //get the neighbors node
            size_t Index1=nIt->nodeID;

            vcg::Point3f pos1;
            pos1.Import<ScalarType>(Graph[Index1]->pos);

            ScalarType CurrGeoL=nIt->geoDistance(max_deviation,penalty_drift);
            assert(MaxGeoL>=MinGeoL);
            vcg::Color4b c=vcg::Color4b::ColorRamp(MaxGeoL,MinGeoL,CurrGeoL);
            vcg::glColor(c);

            vcg::glVertex(pos0);
            vcg::glVertex(pos1);
        }
        glEnd();
        glPopAttrib();
    }

    void GlDrawBarriers()
    {
        std::vector<std::pair<CoordType,CoordType> > barriers;

        for (size_t i=0;i<Graph.size();i++)
        {
            CoordType pos0=Graph[i]->pos;
            for (NeighborsIterator nIt = Graph[i]->neighbors.begin(); nIt != Graph[i]->neighbors.end(); ++nIt)
            {
                if ((*nIt).Valid)continue;
                size_t Index1=nIt->nodeID;
                CoordType pos1=Graph[Index1]->pos;
                barriers.push_back(std::pair<CoordType,CoordType>(pos0,pos1));
            }
        }

        glPushAttrib(GL_ALL_ATTRIB_BITS);
        glDisable(GL_LIGHTING);
        glDepthRange(0,0.99999);

        glLineWidth(2);
        vcg::Color4b c=vcg::Color4b(255,255,0,255);

        glBegin(GL_LINES);
        for (size_t i=0;i<barriers.size();i++)
        {

            CoordType pos0=barriers[i].first;
            CoordType pos1=barriers[i].second;

            vcg::glColor(c);

            vcg::glVertex(pos0);
            vcg::glVertex(pos1);
        }
        glEnd();
        glPopAttrib();
    }

    void GlDrawNodeDebugInfo(size_t Node,
                             ScalarType max_deviation,
                             ScalarType penalty_drift,
                             ScalarType scaleFact=0.0005,
                             bool DrawArrows=true,
                             bool DrawNeigh=true,
                             vcg::Color4b currC=vcg::Color4b(0,255,0,255))
    {
        if (Graph.empty())return;

        //then draw the arrow
        CoordType dir=getDirection(Node);
        vcg::Point3f dirf;
        dirf.Import(dir);
        dirf*=mesh.bbox.Diag()*scaleFact*10;

        vcg::Point3f pos;
        pos.Import<ScalarType>(Graph[Node]->pos);

        //move slightly close to barycenter of face
        if (IsEdgeNode(Node))
        {
            EdgeGeodesicNode* edgeN=static_cast<EdgeGeodesicNode*>(Graph[Node]);
            int faceID=edgeN->faceID;
            vcg::Point3f bary;
            bary.Import((mesh.face[faceID].P(0)+
                         mesh.face[faceID].P(1)+
                         mesh.face[faceID].P(2))/3.0);

            pos=pos*0.9+bary*0.1;
        }

        vcg::glColor(currC);

        float size=mesh.bbox.Diag()*scaleFact;
        vcg::Add_Ons::glPoint<vcg::Add_Ons::DMSolid>(pos,size,5,5);

        if (DrawArrows)
            vcg::Add_Ons::glArrow<vcg::Add_Ons::DMSolid>(pos,pos+dirf,size/2,size*2,size,5,5);

        if (DrawNeigh)
            GlDrawConnections(Node,max_deviation,penalty_drift);
    }

    void GlDrawDisabledLinks(vcg::Color4b currC=vcg::Color4b(255,0,0,255))
    {
        glPushAttrib(GL_ALL_ATTRIB_BITS);
        glDisable(GL_LIGHTING);
        glDepthRange(0,0.99999);
        glLineWidth(2);

        vcg::glColor(currC);

        glBegin(GL_LINES);

        for (size_t Index0=0;Index0<Graph.size();Index0++)
        {
            vcg::Point3f pos0;
            pos0.Import<ScalarType>(Graph[Index0]->pos);
            for (NeighborsIterator nIt = Graph[Index0]->neighbors.begin(); nIt != Graph[Index0]->neighbors.end(); ++nIt)
            {
                if (nIt->Valid)continue;

                //get the neighbors node
                size_t Index1=nIt->nodeID;

                vcg::Point3f pos1;
                pos1.Import<ScalarType>(Graph[Index1]->pos);

                vcg::glVertex(pos0);
                vcg::glVertex(pos1);
            }

        }
        glEnd();
        glPopAttrib();
    }

    void GlDrawSing(ScalarType max_deviation,
                    ScalarType penalty_drift)
    {
        for (size_t i=nEdgeNodes;i<nEdgeNodes+nSingNodes;++i)
        {
            assert(Graph[i]->nodeType==SingularityNodeType);
            GlDrawNodeDebugInfo(i,max_deviation,penalty_drift);
        }
    }

    void GlDrawNodes(bool only_selected=true,
                     bool draw_dir=true)
    {
        glPushAttrib(GL_ALL_ATTRIB_BITS);
        glDisable(GL_LIGHTING);
        glDepthRange(0,0.99999);

        for (size_t i=0;i<Graph.size();++i)
        {
            CoordType pos=Graph[i]->pos;

            if ((only_selected)&&(!Graph[i]->Selected))continue;

            //move slightly close to barycenter of face
            if (IsEdgeNode(i))
            {
                EdgeGeodesicNode* edgeN=static_cast<EdgeGeodesicNode*>(Graph[i]);
                int faceID=edgeN->faceID;
                CoordType bary=(mesh.face[faceID].P(0)+
                                mesh.face[faceID].P(1)+
                                mesh.face[faceID].P(2))/3.0;

                pos=pos*0.9+bary*0.1;
            }

            if (Graph[i]->Selected)
            {
                vcg::glColor(vcg::Color4b(255,0,255,255));
                glPointSize(8);
            }
            else
                if (!Graph[i]->Valid)
                {
                    vcg::glColor(vcg::Color4b(255,0,0,255));
                    glPointSize(4);
                }
                else
                {
                    vcg::glColor(vcg::Color4b(64,64,64,255));
                    glPointSize(4);
                }
            glBegin(GL_POINTS);
            vcg::glVertex(pos);
            glEnd();

            if (!draw_dir)continue;

            CoordType dir=getDirection(i);
            dir.Normalize();
            dir*=mesh.bbox.Diag()*0.002;

            if (Graph[i]->Selected)
            {
                vcg::glColor(vcg::Color4b(255,0,255,255));
                glLineWidth(8);
            }
            else
                if (!Graph[i]->Valid)
                {
                    vcg::glColor(vcg::Color4b(255,0,0,255));
                    glPointSize(4);
                }
                else
                {
                    vcg::glColor(vcg::Color4b(64,64,64,255));
                    glLineWidth(4);
                }
            glBegin(GL_LINES);
            vcg::glVertex(pos);
            vcg::glVertex(pos+dir);
            glEnd();
        }

        glPopAttrib();
    }

    void GlDrawVert(int VIndex)
    {
        glPushAttrib(GL_ALL_ATTRIB_BITS);
        glDisable(GL_LIGHTING);
        glDepthRange(0,0.99999);
        glPointSize(20);
        glBegin(GL_POINTS);
        vcg::glColor(vcg::Color4b(255,255,0,255));
        vcg::glVertex(mesh.vert[VIndex].P());
        glEnd();
        glPopAttrib();
    }
#endif


    bool IsSingularity(size_t IndexN)const
    {
        assert((IndexN>=0)&&(IndexN<Graph.size()));
        return (Graph[IndexN]->nodeType==SingularityNodeType);
    }

    bool IsSink(size_t IndexN)const
    {
        assert((IndexN>=0)&&(IndexN<Graph.size()));
        return (Graph[IndexN]->nodeType==SinkNodeType);
    }

    bool IsEdgeNode(size_t IndexN)const
    {
        assert((IndexN>=0)&&(IndexN<Graph.size()));
        return (Graph[IndexN]->nodeType==EdgeNodeType);
    }

    bool IsBorderNode(size_t IndexN)
    {
        assert((IndexN>=0)&&(IndexN<Graph.size()));
        assert(IsEdgeNode(IndexN));
        EdgeGeodesicNode *edgeN = static_cast<EdgeGeodesicNode*>(Graph[IndexN]);
        return (edgeN->IsBorder);
    }

    int GetVertexSing(size_t IndexN)
    {
        assert(IsSingularity(IndexN)||IsSink(IndexN));
        size_t IndexV=-1;
        if (IsSingularity(IndexN))
        {
            //if it is a singulairity return the relative vertex
            SingularityGeodesicNode* singNode=static_cast<SingularityGeodesicNode*>(Graph[IndexN]);
            IndexV=singNode->vertexID;
        }
        else
        {
            //otherwise go on the other side
            SinkGeodesicNode* sinkNode=static_cast<SinkGeodesicNode*>(Graph[IndexN]);
            size_t OppositeN=sinkNode->oppositeSingularityNode;
            assert(IsSingularity(OppositeN));
            //and then return the relative singularity vertex
            SingularityGeodesicNode* singNode=static_cast<SingularityGeodesicNode*>(Graph[OppositeN]);
            IndexV=singNode->vertexID;
        }
        assert((IndexV>=0)&&(IndexV<mesh.vert.size()));
        assert(IsSingularVert(IndexV));
        return IndexV;
    }

    //void GetFaceDir(const size_t Node0,
    //                size_t &faceID,
    //                size_t &M4Dir)
    //{

    //}


    void GetFaceDir(const size_t Node0,
                    size_t &faceID,
                    size_t &M4Dir)const
    {
        assert(!IsSink(Node0));         //it cannot be a sink
        assert(!IsSingularity(Node0));  //or a singularity

        EdgeGeodesicNode* edgeN=static_cast<EdgeGeodesicNode*>(Graph[Node0]);
        faceID=edgeN->faceID;
        M4Dir=edgeN->M4Dir;
    }

    void GetFaceEdgeDir(const size_t IndexN,
                        size_t &faceID,
                        size_t &edgeIndex,
                        size_t &M4Dir)
    {
        assert(IsEdgeNode(IndexN));
        EdgeGeodesicNode *edgeN = static_cast<EdgeGeodesicNode*>(Graph[IndexN]);
        faceID=edgeN->faceID;
        edgeIndex=edgeN->edgeIdx;
        M4Dir=edgeN->M4Dir;
    }

    void GetFaceDir(const size_t Node0,
                    const size_t Node1,
                    size_t &faceID,
                    size_t &M4Dir)const
    {
        assert(!IsSink(Node0));         //the first one cannot be a sink
        assert(!IsSingularity(Node1));  //the second one cannot be a singularity

        if (IsEdgeNode(Node0))          //in this case is pretty simple
        {
            EdgeGeodesicNode* edgeN=static_cast<EdgeGeodesicNode*>(Graph[Node0]);
            faceID=edgeN->faceID;
            M4Dir=edgeN->M4Dir;
        }
        else
        {
            assert(IsSingularity(Node0));//must be a singularity
            assert(IsEdgeNode(Node1));   //the closest must be an edge
            EdgeGeodesicNode* edgeN=static_cast<EdgeGeodesicNode*>(Graph[Node1]);

            size_t OppNIndex=edgeN->OppositeID;
            assert(IsEdgeNode(OppNIndex));
            EdgeGeodesicNode* OppEdgeN=static_cast<EdgeGeodesicNode*>(Graph[OppNIndex]);
            faceID=OppEdgeN->faceID;
            M4Dir=OppEdgeN->M4Dir;
        }
    }

    void GetNeighInfo(const size_t NodeIndex,
                      std::vector<size_t> &NeighNodes,
                      std::vector<size_t> &FacesId,
                      std::vector<size_t> &M4Dirs,
                      std::vector<bool> &Valid)
    {
        NeighNodes.clear();
        FacesId.clear();
        M4Dirs.clear();
        Valid.clear();
        assert(NodeIndex>=0);
        assert(NodeIndex<Graph.size());

        for (NeighborsIterator NeighIte=Graph[NodeIndex]->neighbors.begin();
             NeighIte!=Graph[NodeIndex]->neighbors.end();NeighIte++)
        {
            size_t FaceId;
            size_t M4Dir;
            size_t Node1=(*NeighIte).nodeID;
            GetFaceDir(NodeIndex,Node1,FaceId,M4Dir);
            NeighNodes.push_back(Node1);
            FacesId.push_back(FaceId);
            M4Dirs.push_back(M4Dir);
            Valid.push_back((*NeighIte).Valid);
        }
    }

    void GetFacesDir(const Path &curr_p,
                     std::vector<std::pair<size_t,size_t> > &FaceDir,
                     bool IsLoop=false)
    {
        FaceDir.clear();
        size_t limitsize=curr_p.nodes.size();
        if (!IsLoop)limitsize--;

        for (size_t i=0;i<limitsize;i++)
        {
            size_t Node0=curr_p.nodes[i];
            size_t Node1=curr_p.nodes[(i+1)%curr_p.nodes.size()];
            assert(Node0!=Node1);
            size_t faceID;
            size_t M4Dir;
            GetFaceDir(Node0,Node1,faceID,M4Dir);
            FaceDir.push_back(std::pair<size_t,size_t>(faceID,M4Dir));
        }
    }

    void GetFacesDirM2(const Path &curr_p,
                       std::vector<std::pair<size_t,size_t> > &FaceDir,
                       bool IsLoop=false)const
    {
        FaceDir.clear();
        size_t limitsize=curr_p.nodes.size();
        if (!IsLoop)limitsize--;

        for (size_t i=0;i<limitsize;i++)
        {
            size_t Node0=curr_p.nodes[i];
            size_t Node1=curr_p.nodes[(i+1)%curr_p.nodes.size()];
            assert(Node0!=Node1);
            size_t faceID;
            size_t M4Dir;
            GetFaceDir(Node0,Node1,faceID,M4Dir);
            FaceDir.push_back(std::pair<size_t,size_t>(faceID,M4Dir));
            M4Dir=(M4Dir+2)%4;
            FaceDir.push_back(std::pair<size_t,size_t>(faceID,M4Dir));
        }
    }

    void GetFacesDirM2Sing(const Path &curr_p,
                           std::vector<std::pair<size_t,size_t> > &FaceDir,
                           bool IsLoop=false)const
    {
        FaceDir.clear();
        size_t limitsize=curr_p.nodes.size();
        if (!IsLoop)limitsize--;

        for (size_t i=0;i<limitsize;i++)
        {
            size_t Node0=curr_p.nodes[i];
            size_t Node1=curr_p.nodes[(i+1)%curr_p.nodes.size()];
            assert(Node0!=Node1);
            size_t faceID;
            size_t M4Dir;
            GetFaceDir(Node0,Node1,faceID,M4Dir);
            FaceDir.push_back(std::pair<size_t,size_t>(faceID,M4Dir%2));
        }
    }

    size_t OppositeSink(int IndexNode)
    {
        assert(IsSink(IndexNode));
        SinkGeodesicNode* sinkN=static_cast<SinkGeodesicNode*>(Graph[IndexNode]);
        return(sinkN->oppositeSingularityNode);
    }

    size_t OppositeSing(int IndexNode)
    {
        assert(IsSingularity(IndexNode));
        SingularityGeodesicNode* singN=static_cast<SingularityGeodesicNode*>(Graph[IndexNode]);
        return(singN->oppositeSinkNode);
    }

    size_t NumNodes()const
    {
        return Graph.size();
    }

    size_t NumFaces()
    {
        return mesh.face.size();
    }

    MeshType& Mesh()const
    {
        return mesh;
    }

    CoordType NodePos(size_t indexN)
    {
        assert(indexN>=0);
        assert(indexN<Graph.size());
        return (Graph[indexN]->pos);
    }

    size_t NodeDir(size_t indexN)
    {
        assert(indexN>=0);
        assert(indexN<Graph.size());
        assert(IsEdgeNode(indexN));
        EdgeGeodesicNode* edgeN = static_cast<EdgeGeodesicNode*>(Graph[indexN]);
        return (edgeN->M4Dir);
    }

    size_t NodeFace(size_t indexN)
    {
        assert(indexN>=0);
        assert(indexN<Graph.size());
        assert(IsEdgeNode(indexN));
        EdgeGeodesicNode* edgeN = static_cast<EdgeGeodesicNode*>(Graph[indexN]);
        return (edgeN->faceID);
    }

    int OppositeNodeSameF(size_t indexN)
    {
        assert(indexN>=0);
        assert(indexN<Graph.size());
        assert(IsEdgeNode(indexN));
        EdgeGeodesicNode* edgeN = static_cast<EdgeGeodesicNode*>(Graph[indexN]);
        return (edgeN->OppositeSameFID);
    }

    int OppositeNode(size_t indexN)
    {
        assert(indexN>=0);
        assert(indexN<Graph.size());
        assert(IsEdgeNode(indexN));
        EdgeGeodesicNode* edgeN = static_cast<EdgeGeodesicNode*>(Graph[indexN]);
        return (edgeN->OppositeID);
    }

    void OppositeNode(const std::vector<size_t> &indexN,
                      std::vector<size_t> &indexOpp)
    {
        for (size_t i=0;i<indexN.size();i++)
            indexOpp.push_back(OppositeNode(indexN[i]));
    }

    void InvertPath(Path &To_Invert)
    {
        //std::cout << "size ... " << To_Invert.nodes.size() << std::endl;
        for (size_t i=0;i<To_Invert.nodes.size();i++)
        {
            //std::cout << "i ... " << i << std::endl;

            size_t currN=To_Invert.nodes[i];
            assert(currN>=0);
            assert(currN<Graph.size());
            //the first one can be a border one in this case no need to invert
            if ((i==0) && IsEdgeNode(currN) && (IsBorderNode(currN)))To_Invert.nodes[i]=currN;
            else
                //otherwise invert
                if (IsEdgeNode(currN))To_Invert.nodes[i]=OppositeNode(currN);
                else
                    if (IsSink(currN))To_Invert.nodes[i]=OppositeSink(currN);
                    else
                        if (IsSingularity(currN))To_Invert.nodes[i]=OppositeSing(currN);
            //final check
            assert(To_Invert.nodes[i]!=-1);
            //        assert(To_Invert.nodes[i]>=0);
            //        assert(To_Invert.nodes[i]<Graph.size());
        }
        std::reverse(To_Invert.nodes.begin(),To_Invert.nodes.end());
    }

    void TangentialNodes(size_t indexN,
                         std::vector<size_t> &Tangential,
                         int N_dir=2)
    {
        assert(indexN>=0);
        assert(indexN<Graph.size());
        assert(IsEdgeNode(indexN));
        Tangential.clear();
        EdgeGeodesicNode* edgeN = static_cast<EdgeGeodesicNode*>(Graph[indexN]);
        CoordType pos=edgeN->pos;
        size_t testF=edgeN->faceID;
        NodeListConstIterator First,Last;
        GetCoordNodes(pos,First,Last);
        for (NodeListConstIterator Curr=First;Curr!=Last;Curr++)
        {
            int CurrNodeI=(*Curr);
            assert(IsEdgeNode(CurrNodeI));
            EdgeGeodesicNode* Curr_edgeN = static_cast<EdgeGeodesicNode*>(Graph[CurrNodeI]);

            size_t TestM4Dir=edgeN->M4Dir;
            if (Curr_edgeN->faceID!=testF)
                TestM4Dir=vcg::tri::CrossField<MeshType>::FollowDirection(mesh.face[testF],mesh.face[Curr_edgeN->faceID],TestM4Dir);

            TestM4Dir=TestM4Dir%N_dir;
            if ((Curr_edgeN->M4Dir % N_dir)!=(int)TestM4Dir)Tangential.push_back(CurrNodeI);
        }
    }

    template <class MESH_TYPE>
    class SelEdgePred
    {
    public:
        bool operator()(vcg::face::Pos<typename MESH_TYPE::FaceType> ep)
        {
            VertexType *v0=ep.f->V0(ep.z);
            VertexType *v1=ep.f->V1(ep.z);
            if (v0->IsV()&&v1->IsV())return true;
            if (v0->IsB()&&v1->IsV())return true;
            if (v0->IsV()&&v1->IsB())return true;
            return false;
            //return (v0->IsV()&&v1->IsV());
        }
    };


    template<class MESH_TYPE>
    class SelMidPointFunctor : public std::unary_function<vcg::face::Pos<typename MESH_TYPE::FaceType> , typename MESH_TYPE::CoordType>
    {
    public:

        void operator()(typename MESH_TYPE::VertexType &nv, const vcg::face::Pos<typename MESH_TYPE::FaceType> &ep){
            CoordType p0=ep.f->V0(ep.z)->P();
            CoordType p1=ep.f->V1(ep.z)->P();
            nv.P()=p1*0.5 + p0*0.5;
        }

        template<class FL_TYPE>
        vcg::TexCoord2<FL_TYPE,1> WedgeInterp(vcg::TexCoord2<FL_TYPE,1> &t0,
                                              vcg::TexCoord2<FL_TYPE,1> &t1)
        {
            vcg::TexCoord2<FL_TYPE,1> tmp;
            assert(t0.n()== t1.n());
            tmp.n()=t0.n();
            tmp.t()=(t0.t()+t1.t())/2.0;
            return tmp;
        }

    };

    void GetCoordNodes(const CoordType &Pos,
                       NodeListConstIterator &First,
                       NodeListConstIterator &End)
    {
        assert(CoordNodes.count(Pos)>0);
        First=CoordNodes[Pos].begin();
        End=CoordNodes[Pos].end();
    }


    VertexType* GetSingVertex(size_t &NodeI)
    {
        assert((NodeI>=0)&&(NodeI<Graph.size()));
        assert(IsSingularity(NodeI)||IsSink(NodeI));
        if (IsSingularity(NodeI))
        {
            SingularityGeodesicNode* SingN=static_cast<SingularityGeodesicNode*>(Graph[NodeI]);
            size_t singVIndex=SingN->vertexID;
            return (&mesh.vert[singVIndex]);
        }
        else
        {
            assert(IsSink(NodeI));
            SinkGeodesicNode* SinkN=static_cast<SinkGeodesicNode*>(Graph[NodeI]);
            int OppositeS=SinkN->oppositeSingularityNode;
            assert(OppositeS>=0);
            assert(OppositeS<NumNodes());
            SingularityGeodesicNode* SingN=static_cast<SingularityGeodesicNode*>(Graph[OppositeS]);
            size_t singVIndex=SingN->vertexID;
            return (&mesh.vert[singVIndex]);
        }
    }

    //return all possible direction for a given pos
    //considering %2.. need to check if better one face or the other
    //it is not complete
    void GetPossibleDirNodes(const CoordType &Pos,
                             std::vector<size_t> &dirNodes,
                             int NumDir=2)
    {
        dirNodes.clear();

        assert(CoordNodes.count(Pos)>0);

        NodeListConstIterator First,End;
        GetCoordNodes(Pos,First,End);

        //if it is a singularity then add all of them
        if (IsSingularity(*First)||(IsSink(*First)))
        {
            dirNodes=std::vector<size_t>(First,End);
            return;
        }
        //then get one per possible direction per Face
        std::vector<bool> FoundDir(NumDir,false);


        //get the first face it belongs to
        size_t FirstNodeI=(*First);
        EdgeGeodesicNode* EdgeN=static_cast<EdgeGeodesicNode*>(Graph[FirstNodeI]);
        size_t IndexF=EdgeN->faceID;
        for (NodeListConstIterator Curr=First;Curr!=End;Curr++)
        {
            size_t NodeI=(*Curr);
            EdgeGeodesicNode* currN=static_cast<EdgeGeodesicNode*>(Graph[NodeI]);
            //if (currN->faceID!=IndexF)continue;
            int M4DirTest=currN->M4Dir;
            int NextF=currN->faceID;
            CoordType DirNode=getDirection(NodeI);
            //transform back to the original face
            if (NextF!=(int)IndexF)
                M4DirTest=vcg::tri::CrossField<MeshType>::FollowDirectionI(mesh.face[NextF],
                                                                           mesh.face[IndexF],
                                                                           DirNode);

            size_t currD=(M4DirTest % NumDir);

            //insert only one direction %Ndir
            if (FoundDir[currD])continue;
            dirNodes.push_back(NodeI);
            FoundDir[currD]=true;
        }
    }

    //    void GetEdgeNodeFromPosDir(const CoordType &Pos,CoordType &Dir,std::vector<size_t> &Nodes)
    //    {
    //        NodeListConstIterator Curr,First,End;
    //        GetCoordNodes(Vpos,First,End);
    //        for (Curr=First;First!=End;First++)
    //        {
    //            if ((*Curr).NodeType!=EdgeNodeType)continue;
    //            EdgeGeodesicNode* CurrN=static_cast<EdgeGeodesicNode*>(Graph[NodeI]);
    //            size_t IndexN= CurrN->faceID;
    //            int     edgeIdx;        ///< The index of the edge the node belongs to.
    //            int     M4Dir;          ///< The number of the vector of the cross field of M4.
    //            int     OppositeID;     ///< The index of the opposite Node.
    //            int     OppositeSameFID;///< The index of the opposite Node  considering the same Face
    //            bool    IsBorder;
    //        }
    //    }

    void InterpolationFace(size_t IndexN,
                           size_t &FaceIndex,
                           CoordType &bary)
    {
        CoordType pos;
        if (IsEdgeNode(IndexN))          //in this case is pretty simple
        {
            EdgeGeodesicNode* edgeN=static_cast<EdgeGeodesicNode*>(Graph[IndexN]);
            FaceIndex=edgeN->faceID;
            pos=edgeN->pos;
        }
        else
            if (IsSingularity(IndexN))
            {
                SingularityGeodesicNode *singN=static_cast<SingularityGeodesicNode*>(Graph[IndexN]);
                FaceIndex=singN->exitFace;
                pos=singN->pos;
            }
            else
            {
                assert(IsSink(IndexN));
                int OppN=OppositeSink(IndexN);
                assert(IsSingularity(OppN));
                SingularityGeodesicNode* singN=static_cast<SingularityGeodesicNode*>(Graph[OppN]);
                FaceIndex=singN->exitFace;
                pos=singN->pos;
            }
        assert(FaceIndex>=0);
        assert(FaceIndex<mesh.face.size());
        FaceType *f=&mesh.face[FaceIndex];
        vcg::InterpolationParameters(*f,pos,bary);
    }

private:

    bool SplitVEdge()
    {
        bool hasSingular = vcg::tri::HasPerVertexAttribute(mesh,std::string("Singular"));
        bool hasSingularIndex = vcg::tri::HasPerVertexAttribute(mesh,std::string("SingularIndex"));

        assert(hasSingular);
        assert(hasSingularIndex);
        typename MeshType::template PerVertexAttributeHandle<bool> Handle_Singular;
        typename MeshType::template PerVertexAttributeHandle<int> Handle_SingularIndex;

        Handle_Singular=vcg::tri::Allocator<MeshType>::template GetPerVertexAttribute<bool>(mesh,std::string("Singular"));
        Handle_SingularIndex=vcg::tri::Allocator<MeshType>::template GetPerVertexAttribute<int>(mesh,std::string("SingularIndex"));

        //save all singularities and valences
        std::map<CoordType,size_t> SingValences;
        for (size_t i=0;i<mesh.vert.size();i++)
        {
            if (!Handle_Singular[i])continue;
            int SingIndex=Handle_SingularIndex[i];
            CoordType CurrPos=mesh.vert[i].P();
            SingValences[CurrPos]=SingIndex;
        }

        //then Split the edges
        SelEdgePred<MeshType> SelEP;
        SelMidPointFunctor<MeshType> MidEP;
        bool refined=vcg::tri::RefineE<MeshType,SelMidPointFunctor<MeshType>,SelEdgePred<MeshType> >(mesh,MidEP,SelEP);
        vcg::tri::UpdateTopology<MeshType>::FaceFace(mesh);
        vcg::tri::UpdateTopology<MeshType>::VertexFace(mesh);

        //then set singularities again
        for (size_t i=0;i<mesh.vert.size();i++)
        {
            CoordType CurrPos=mesh.vert[i].P();
            if (SingValences.count(CurrPos)==0)
            {
                Handle_Singular[i]=false;
                Handle_SingularIndex[i]=0;
            }else
            {
                Handle_Singular[i]=true;
                int SingIndex=SingValences[CurrPos];
                Handle_SingularIndex[i]=SingIndex;
            }
        }

        //vcg::tri::CrossField<MeshType>::UpdateSingularByCross(mesh);
        //mesh.InitSingFlag();

        //update normals
        vcg::tri::UpdateNormal<MeshType>::PerFaceNormalized(mesh);
        vcg::tri::UpdateNormal<MeshType>::PerVertexNormalized(mesh);
        return refined;
    }

    bool SplitSing(size_t multCard=2)
    {
        std::vector<size_t> Card(mesh.vert.size(),0);
        //update the cardinality of each singularity (number of faces)
        for (size_t i=0;i<mesh.face.size();i++)
        {
            for (size_t j=0;j<3;j++)
            {
                int IndexV=vcg::tri::Index(mesh,mesh.face[i].V(j));
                Card[IndexV]++;
            }
        }

        //select all edges that have at least one singularity
        vcg::tri::UpdateFlags<MeshType>::VertexClearV(mesh);
        for (size_t i=0;i<mesh.face.size();i++)
        {
            for (size_t j=0;j<3;j++)
            {
                int IndexV0=vcg::tri::Index(mesh,mesh.face[i].V0(j));
                int IndexV1=vcg::tri::Index(mesh,mesh.face[i].V1(j));
                //not possible two nearby singularities
                assert(!(IsSingularVert(IndexV0)&&IsSingularVert(IndexV1)));
                int to_test=-1;
                int to_split=-1;

                //check singularity and valence
                if (IsSingularVert(IndexV0)&&(!IsSingularVert(IndexV1)))
                {
                    to_test=IndexV0;
                    to_split=IndexV1;
                }
                if (IsSingularVert(IndexV1)&&(!IsSingularVert(IndexV0)))
                {
                    to_test=IndexV1;
                    to_split=IndexV0;
                }
                if (to_test==-1)continue;

                int ExpVal=expectedValence(to_test);
                int NumF=Card[to_test];

                if (ExpVal<(NumF*multCard))
                    mesh.vert[to_split].SetV();
            }
        }
        size_t Faces0=mesh.fn;
        //then Split the edges
        bool refined=SplitVEdge();
        size_t Faces1=mesh.fn;
        std::cout << "Added Faces " << Faces1 - Faces0 << std::endl;

        return refined;
    }

    void SplitAdjSing()
    {
        vcg::tri::UpdateFlags<MeshType>::VertexClearV(mesh);
        for (size_t i=0;i<mesh.vert.size();i++)
            if (IsSingularVert(i))mesh.vert[i].SetV();

        //then Split the edges
        SplitVEdge();
    }

public:

    void  WeightLenghtByVertexQ()
    {
        //        assert((perc>=0)&&(perc<=1));

        //        //first find percentily of quality
        //        std::vector<ScalarType > QVect;
        //        for (size_t i=0;i<mesh.vert.size();i++)
        //            QVect.push_back(mesh.vert[i].Q());

        //        //then find percentiles
        //        std::sort(QVect.begin(),QVect.end());
        //        int Index0=floor((ScalarType)QVect.size()*perc+0.5);
        //        int Index1=floor((ScalarType)QVect.size()*(1-perc)+0.5);

        //        //then find limits
        //        ScalarType minQ=QVect[Index0];
        //        ScalarType maxQ=QVect[Index1];

        //then set lenght
        for (size_t i=0;i<Graph.size();i++)
        {
            size_t FaceIndex0;
            CoordType bary0;
            InterpolationFace(i,FaceIndex0,bary0);
            FaceType *f0=&mesh.face[FaceIndex0];

            //then interpolate quality
            ScalarType Q0=f0->V(0)->Q()*bary0.X()+
                    f0->V(1)->Q()*bary0.Y()+
                    f0->V(2)->Q()*bary0.Z();

            for (NeighborsIterator NeighIte=Graph[i]->neighbors.begin();
                 NeighIte!=Graph[i]->neighbors.end();NeighIte++)
            {
                size_t IndexN=(*NeighIte).nodeID;
                size_t FaceIndex1;
                CoordType bary1;
                InterpolationFace(IndexN,FaceIndex1,bary1);
                FaceType *f1=&mesh.face[FaceIndex1];

                //then interpolate quality
                ScalarType Q1=f1->V(0)->Q()*bary1.X()+
                        f1->V(1)->Q()*bary1.Y()+
                        f1->V(2)->Q()*bary1.Z();

                ScalarType AvgQ=(Q0+Q1)/2;
                //then crop It
                //                if (AvgQ<minQ)AvgQ=minQ;
                //                if (AvgQ>maxQ)AvgQ=maxQ;
                //                //and set it
                (*NeighIte).Weight=1/AvgQ;
            }
        }
    }

    //    void WeightByCurvatureAnisotropy()
    //    {
    //        FieldSmoother<MeshType>::InitByCurvature(mesh,3,true);
    //        WeightLenghtByVertexQ();
    //        vcg::tri::UpdateColor<MeshType>::PerVertexQualityRamp(mesh);
    //        vcg::tri::UpdateColor<MeshType>::PerFaceFromVertex(mesh);
    //    }

    void WeightUniform()
    {
        for (int i=0;i<Graph.size();i++)
        {
            for (NeighborsIterator NeighIte=Graph[i]->neighbors.begin();
                 NeighIte!=Graph[i]->neighbors.end();NeighIte++)
                (*NeighIte).Weight=1;
        }
        vcg::tri::UpdateColor<MeshType>::PerFaceConstant(mesh,vcg::Color4b(200,200,200,255));
    }

    void PathPos(const Path &Path,
                 std::vector<CoordType> &PathPos)
    {
        PathPos.clear();

        for (size_t i=0;i<Path.nodes.size();++i)
            PathPos.push_back(Graph[Path.nodes[i]]->pos);

    }

    void SetValid(size_t IndexN,bool IsValid)
    {
        Graph[IndexN]->Valid=IsValid;
    }

    bool IsValid(size_t IndexN)
    {
        return (Graph[IndexN]->Valid);
    }

    void SetNeighValidity(size_t IndexN,
                          const std::vector<bool> &IsValid)
    {
        size_t curr_i=0;
        for (NeighborsIterator NeighIte=Graph[IndexN]->neighbors.begin();
             NeighIte!=Graph[IndexN]->neighbors.end();NeighIte++)
        {
            assert(curr_i<IsValid.size());
            (*NeighIte).Valid=IsValid[curr_i];
            curr_i++;
        }
    }

    void InvalidateArcsOnNonValidNodes()
    {
        //first set as unvalid exiting arcs
        for (size_t i=0;i<Graph.size();i++)
        {
            if (Graph[i]->Valid)continue;
            for (NeighborsIterator NeighIte=Graph[i]->neighbors.begin();
                 NeighIte!=Graph[i]->neighbors.end();NeighIte++)
                (*NeighIte).Valid=false;
        }

        //then set as unvalid entering arcs
        for (size_t i=0;i<Graph.size();i++)
        {
            if (!Graph[i]->Valid)continue;
            for (NeighborsIterator NeighIte=Graph[i]->neighbors.begin();
                 NeighIte!=Graph[i]->neighbors.end();NeighIte++)
            {
                size_t IndexN=(*NeighIte).nodeID;
                if (Graph[IndexN]->Valid)continue;
                (*NeighIte).Valid=false;
            }
        }

    }

    void SetAllValid()
    {
        //first set as unvalid exiting arcs
        for (int i=0;i<(int)Graph.size();i++)
        {
            Graph[i]->Valid=true;
            for (NeighborsIterator NeighIte=Graph[i]->neighbors.begin();
                 NeighIte!=Graph[i]->neighbors.end();NeighIte++)
                (*NeighIte).Valid=true;
        }
    }

    void SetValidSink(const std::vector<size_t> &SinkIndex)
    {
        //first set as unvalid all sink
        for (int i=0;i<Graph.size();i++)
        {
            if (!IsSink(i))continue;
            Graph[i]->Valid=false;
        }
        for (int i=0;i<SinkIndex.size();i++)
        {
            size_t IndexN=SinkIndex[i];
            assert(IsSink(IndexN));
            Graph[i]->Valid=true;
        }
        InvalidateArcsOnNonValidNodes();
    }


    void SetAllInvalid()
    {
        //first set as unvalid exiting arcs
        for (int i=0;i<Graph.size();i++)
        {
            Graph[i]->Valid=false;
            //            for (NeighborsIterator NeighIte=Graph[i]->neighbors.begin();
            //                 NeighIte!=Graph[i]->neighbors.end();NeighIte++)
            //                (*NeighIte).Valid=false;
        }
    }

    ScalarType TotalWeight(ScalarType max_deviation,
                           ScalarType penalty_drift)
    {
        ScalarType TotW=0;
        //first set as unvalid exiting arcs
        for (int i=0;i<Graph.size();i++)
        {
            Graph[i]->Valid=true;
            for (NeighborsIterator NeighIte=Graph[i]->neighbors.begin();
                 NeighIte!=Graph[i]->neighbors.end();NeighIte++)
                TotW+=NeighIte->geoDistance(max_deviation, penalty_drift);
        }
        return TotW;
    }

    void SelectNode(size_t IndexN)
    {
        Graph[IndexN]->Selected=true;
    }

    void DeselectAllNodes()
    {
        //first set as unvalid exiting arcs
        for (int i=0;i<Graph.size();i++)
            Graph[i]->Selected=false;
    }

    void GetClosest(const CoordType &TestPos,
                    std::vector<size_t> &IndexNodes,
                    bool excludesinks=true)
    {
        if (Graph.size()==0)return;
        ScalarType minD=std::numeric_limits<ScalarType>::max();
        CoordType minPos;
        for (int i=0;i<Graph.size();i++)
        {
            ScalarType TestDist=(TestPos-Graph[i]->pos).Norm();
            if (TestDist>=minD)continue;
            minD=TestDist;
            minPos=Graph[i]->pos;
        }
        IndexNodes.clear();
        NodeListConstIterator First,End;
        GetCoordNodes(minPos,First,End);
        if (excludesinks)
        {
            for (NodeListConstIterator Curr=First;Curr!=End;Curr++)
            {
                if (IsSink(*Curr))continue;
                IndexNodes.push_back(*Curr);
            }
        }
        else
            IndexNodes.insert(IndexNodes.end(),First,End);
    }

    void TangentialNodes(const std::vector<Path> &Paths,
                         std::vector<size_t> &TangNodes,
                         size_t Ndir)
    {
        TangNodes.clear();
        //then get all the elements from the paths and get the
        //tangential node
        for (size_t i=0;i<Paths.size();i++)
            for (size_t j=0;j<Paths[i].nodes.size();j++)
            {
                size_t currN=Paths[i].nodes[j];
                if(IsSingularity(currN))continue;
                if(IsSink(currN))continue;
                //find the tangential directions
                std::vector<size_t> TangentialN;
                TangentialNodes(currN,TangentialN,Ndir);
                for (size_t k=0;k<TangentialN.size();k++)
                {
                    //may be not valid
                    int TangN=TangentialN[k];
                    if (!IsValid(TangN))continue;
                    SelectNode(TangN);
                    //then insert in the possible sources
                    TangNodes.push_back(TangN);
                }
            }
    }

    void GetBorderNodes(std::vector<size_t> &BorderNodes)
    {
        BorderNodes.clear();
        for (int i=0;i<Graph.size();i++)
        {
            if (!IsValid(i))continue;
            if (!IsEdgeNode(i))continue;
            if (!IsBorderNode(i))continue;
            BorderNodes.push_back(i);
        }
    }

    //    void WeightByTangentDist(std::vector<Path> &Paths,
    //                             ScalarType max_deviation,
    //                             ScalarType penalty_drift,
    //                             size_t Ndir)
    //    {
    //        std::vector<size_t> TangNodes;
    //        TangentialNodes(Paths,TangNodes,Ndir);
    //        DijkstraParam Param;
    //        Param.maxDist=std::numeric_limits<ScalarType>::max();
    //        Param.maxSink=std::numeric_limits<size_t>::max();
    //        Param.sourceGroups.resize(1,TangNodes);

    //        computePerNodeDijsktra(Param,max_deviation,penalty_drift);
    //    }

    void GetPosPairs(std::vector<Path> &Paths,
                     std::vector<bool> &LoopPaths,
                     std::vector<std::vector<std::pair<CoordType,CoordType> > > &EdgePos)
    {
        EdgePos.clear();
        EdgePos.resize(Paths.size());
        for (size_t i=0;i<Paths.size();i++)
        {
            std::vector<CoordType> PerPathPos;

            bool IsLoop=LoopPaths[i];

            //then get the positions
            PathPos(Paths[i],PerPathPos);

            //then add all the paths
            int limit=PerPathPos.size()-1;
            if (IsLoop)
                limit++;

            for (int j=0;j<limit;j++)
            {
                //initialize the edge
                CoordType pos0=PerPathPos[j];
                CoordType pos1=PerPathPos[(j+1)%PerPathPos.size()];

                std::pair<CoordType,CoordType> key=std::pair<CoordType,CoordType>(std::min(pos0,pos1),std::max(pos0,pos1));

                EdgePos[i].push_back(key);
            }
        }
    }

    template <class EdgeMesh>
    void GetEdgeMesh(const std::vector<Path> &splitPaths,
                     const std::vector<bool> &LoopPaths,
                     EdgeMesh &EdgeM,
                     MeshType &tri_mesh)
    {
        EdgeM.Clear();

        //get all the edges
        std::vector<std::pair<CoordType,CoordType> >  EdgePos;
        //this is for faces and directions
        std::vector<std::pair<size_t,size_t> > FaceDir;
        //and the path for each edge
        std::vector<int>  PathIndex;

        assert(LoopPaths.size()==splitPaths.size());
        //for each path create the paths
        for (size_t i=0;i<splitPaths.size();i++)
        {
            std::vector<CoordType> PerPathPos;

            //get faces + directions
            std::vector<std::pair<size_t,size_t> > FaceDirTemp;
            bool IsLoop=LoopPaths[i];
            GetFacesDir(splitPaths[i],FaceDirTemp,IsLoop);

            //then get the positions
            PathPos(splitPaths[i],PerPathPos);

            //safety check
            if (!IsLoop)
            {
                assert(FaceDirTemp.size()==PerPathPos.size()-1);
            }
            else
            {
                assert(FaceDirTemp.size()==PerPathPos.size());
            }
            //then add all the paths
            int limit=PerPathPos.size()-1;
            if (IsLoop)
                limit++;

            for (int j=0;j<limit;j++)
            {
                //initialize the edge
                CoordType pos0=PerPathPos[j];
                CoordType pos1=PerPathPos[(j+1)%PerPathPos.size()];
                if (pos0==pos1)continue;
                std::pair<CoordType,CoordType> key=std::pair<CoordType,CoordType>(std::min(pos0,pos1),std::max(pos0,pos1));

                EdgePos.push_back(key);

                //put the path that the separatrix comes from
                PathIndex.push_back(i);

                //add other info for the split wich will follows
                FaceDir.push_back(FaceDirTemp[j]);
            }
            if (IsLoop)
            {
                assert(FaceDirTemp.size()==PerPathPos.size());
            }
            else
            {
                assert(FaceDirTemp.size()==(PerPathPos.size()-1));
            }
        }

        //create the edge mesh
        EdgeM.InitFromCoords(EdgePos,FaceDir,PathIndex);
        //        std::cout<<"size "<<EdgeM.edge.size()<<std::endl;
        //        std::cout<<"size 0 "<<EdgePos.size()<<std::endl;
        //        std::sort(EdgePos.begin(),EdgePos.end());
        //        auto last = std::unique(EdgePos.begin(), EdgePos.end());
        //        EdgePos.erase(last, EdgePos.end());
        //        std::cout<<"size 1 "<<EdgePos.size()<<std::endl;
        //        assert(0);

        //get all singularities position
        std::set<CoordType> SingPos;
        for (size_t i=0;i<tri_mesh.vert.size();i++)
        {
            //if (!tri_mesh.IsSingularVert(tri_mesh.vert[i]))continue;
            if (!IsSingularVert(i))continue;
            //if (vcg::tri::CrossField<MeshType>::IsSingular(tri_mesh,tri_mesh.vert[i]))continue;
            SingPos.insert(tri_mesh.vert[i].P());
        }
        //set singular vert
        for (size_t i=0;i<EdgeM.vert.size();i++)
        {
            if (SingPos.count(EdgeM.vert[i].P())==0)continue;
            EdgeM.vert[i].IsSing=true;
        }
        //then create the map of faces to edges
        EdgeM.SplitSelfIntersections(tri_mesh);

        EdgeM.InitColorByPath();
    }

    void GetSingExtremes(Path &CurrPath,
                         size_t &SingNode0,
                         size_t &SingNode1)
    {
        //get the two singularities
        //which are the first and the last of the path
        SingNode0=CurrPath.nodes.front();
        size_t SinkNode1=CurrPath.nodes.back();

        //check the first must be a singulatity
        assert(IsSingularity(SingNode0));
        //check the last one must be a sink
        assert(IsSink(SinkNode1));

        //only insert in one sense the index of the nodes MUST be different
        assert(SingNode0!=SinkNode1);
        SingNode1=OppositeSink(SinkNode1);
        assert(IsSingularity(SingNode1));
    }
};
#endif //ANISOTROPIC_GEODESIC_GRAPH
