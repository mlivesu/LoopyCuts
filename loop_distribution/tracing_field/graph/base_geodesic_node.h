#ifndef BASE_NODE_ANISOTROPIC_GEODESIC
#define BASE_NODE_ANISOTROPIC_GEODESIC

#include "neigh_info.h"
#include <list>
#include <deque>

    /**
     * @brief The NodeType enum defines the types of the nodes.
     */
    enum NodeType {
        EdgeNodeType,
        SingularityNodeType,
        SinkNodeType
    };

    /**
     * @brief The BaseGeodesicNode class
     */
    template <class CoordType>
    class BaseGeodesicNode
    {
        typedef typename CoordType::ScalarType ScalarType;

    public:

        typedef std::list<NeighborInfo<ScalarType> >    NeighborsList;
        typedef typename NeighborsList::iterator        NeighborsIterator;
        typedef typename NeighborsList::const_iterator  NeighborsConstIterator;


        NodeType                    nodeType;   ///< The type of node.
        CoordType                   pos;        ///< The current position.
        NeighborsList               neighbors;  ///< The list of neighbors.
        bool                        Selected;
        bool                        Valid;

        //ScalarType NodeW;

        /**
         * @brief BaseGeodesicNode Constructor.
         * @param _nodeType
         * @param _pos
         */
        BaseGeodesicNode(const NodeType _nodeType, const CoordType &_pos)
            : nodeType(_nodeType), pos(_pos) {Selected=false;Valid=true;/*NodeW=1;*/ }

        /**
         * @brief BaseGeodesicNode Copy constructor.
         * @param bgn
         */
        BaseGeodesicNode(const BaseGeodesicNode &bgn)
            : nodeType(bgn.nodeType),
              pos(bgn.pos),
              neighbors(bgn.neighbors),
              Selected(bgn.Selected),
              Valid(bgn.Valid)/*,
                  NodeW(bgn.NodeW)*/{ }

        /**
         * @brief operator = Copy operator.
         * @param bgn
         * @return
         */
        BaseGeodesicNode & operator =(const BaseGeodesicNode &bgn) {
            if (this != &bgn) {
                nodeType = bgn.nodeType;
                pos = bgn.pos;
                neighbors = bgn.neighbors;
                Selected = bgn.Selected;
                Valid = bgn.Valid;
                //NodeW= bgn.NodeW;
            }
            return *this;
        }

#if __cplusplus >= 201103L
        /**
         * @brief BaseGeodesicNode Move constructor.
         * @param bgn
         */
        BaseGeodesicNode(BaseGeodesicNode &&bgn)
            : nodeType(bgn.nodeType),
              pos(bgn.pos),
              neighbors(std::move(bgn.neighbors)),
              Selected(bgn.Selected),
              Valid(bgn.Valid)/*,
                  NodeW(bgn.NodeW)*/{ }

        /**
          * @brief operator = Move operator.
          * @param bgn
          * @return
          */
        BaseGeodesicNode & operator =(BaseGeodesicNode &&bgn)
        {
            if (this != &bgn)
            {
                nodeType = bgn.nodeType;
                pos = bgn.pos;
                neighbors = std::move(bgn.neighbors);
                Selected = bgn.Selected;
                Valid=bgn.Valid;
                //NodeW=bgn.NodeW;
            }
            return *this;
        }
#endif

        /**
         * @brief MinMaxDistance
         * @param minD
         * @param maxD
         * @param max_angle
         * @param drift_penalty
         */
        void MinMaxDistance(ScalarType &minD,
                            ScalarType &maxD,
                            const ScalarType max_angle = 45,
                            const ScalarType drift_penalty = 2)
        {
            maxD = ScalarType(0);
            minD = std::numeric_limits<ScalarType>::max();
            for (NeighborsIterator it = neighbors.begin(); it != neighbors.end(); ++it) {
                ScalarType currD = it->geoDistance(max_angle, drift_penalty);
                if (currD > maxD)
                    maxD = currD;
                if (currD < minD)
                    minD = currD;
            }
        }

        void MinMaxAngle(ScalarType &minA,
                         ScalarType &maxA) {
            maxA = ScalarType(0);
            minA = std::numeric_limits<ScalarType>::max();
            for (NeighborsIterator it = neighbors.begin(); it != neighbors.end(); ++it) {
                if (it->angle > maxA)
                    maxA = it->angle;
                if (it->angle < minA)
                    minA = it->angle;
            }
        }
    };

#endif //BASE_NODE_ANISOTROPIC_GEODESIC
