#ifndef SINK_NODE_ANISOTROPIC_GEODESIC
#define SINK_NODE_ANISOTROPIC_GEODESIC

#include "base_geodesic_node.h"

    /**
     * @brief The SinkGeodesicNode class
     */
    template <class CoordType>
    class SinkGeodesicNode : public BaseGeodesicNode<CoordType>
    {
        typedef BaseGeodesicNode<CoordType> BaseGeodesicNode;

    public:
        size_t  oppositeSingularityNode;    ///< The index of the opposite singularoty node of the sink.

        /**
         * @brief SinkGeodesicNode Constructor.
         * @param _oppositeSingularityNode
         * @param _pos
         */
        SinkGeodesicNode(const size_t _oppositeSingularityNode,
                         const CoordType &_pos)
            : BaseGeodesicNode(SinkNodeType, _pos),
              oppositeSingularityNode(_oppositeSingularityNode) { }

        /**
         * @brief SinkGeodesicNode Copy constructor.
         * @param sgn
         */
        SinkGeodesicNode(const SinkGeodesicNode &sgn)
            : BaseGeodesicNode(sgn),
              oppositeSingularityNode(sgn.oppositeSingularityNode) { }

        /**
         * @brief operator = Copy operator.
         * @param sgn
         * @return
         */
        SinkGeodesicNode & operator =(const SinkGeodesicNode &sgn) {
            if (this != &sgn) {
                BaseGeodesicNode::operator =(sgn);
                oppositeSingularityNode = sgn.oppositeSingularityNode;
            }
            return *this;
        }

#if __cplusplus >= 201103L
        /**
         * @brief SinkGeodesicNode Move constructor.
         * @param sgn
         */
        SinkGeodesicNode(SinkGeodesicNode &&sgn)
            : BaseGeodesicNode(std::move(sgn)),
              oppositeSingularityNode(sgn.oppositeSingularityNode) { }

        /**
         * @brief operator = Move operator.
         * @param sgn
         * @return
         */
        SinkGeodesicNode & operator =(SinkGeodesicNode &&sgn) {
            if (this != &sgn) {
                BaseGeodesicNode::operator =(std::move(sgn));
                oppositeSingularityNode = sgn.oppositeSingularityNode;
            }
            return *this;
        }
#endif
    };
#endif //SINK_NODE_ANISOTROPIC_GEODESIC
