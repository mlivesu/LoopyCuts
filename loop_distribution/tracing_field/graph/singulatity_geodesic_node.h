#ifndef SINGULARITY_NODE_ANISOTROPIC_GEODESIC
#define SINGULARITY_NODE_ANISOTROPIC_GEODESIC

#include "base_geodesic_node.h"

/**
     * @brief The SingularityGeodesicNode class
     */
    template <class CoordType>
    class SingularityGeodesicNode : public BaseGeodesicNode<CoordType>
    {
    public:
        typedef BaseGeodesicNode<CoordType> BaseGeodesicNode;

        size_t      vertexID;           ///the index of the vertex the node belongs to.
        CoordType   mainDirection;      ///the main direction of the separatrix.
        size_t      exitFace;           ///the face the separatrix exit from
        size_t      oppositeSinkNode;   ///the index of the opposite sink node

        /**
         * @brief SingularityGeodesicNode Constructor.
         * @param _vertexID
         * @param _pos
         * @param _mainDirection
         */
        SingularityGeodesicNode(const size_t _vertexID,
                                const CoordType &_pos,
                                const CoordType &_mainDirection,
                                const size_t &_exitFace)
            : BaseGeodesicNode(SingularityNodeType, _pos),
              vertexID(_vertexID),
              mainDirection(_mainDirection),
              exitFace(_exitFace)
        { oppositeSinkNode=std::numeric_limits<size_t>::max();}

        /**
         * @brief SingularityGeodesicNode Copy constructor.
         * @param vgn
         */
        SingularityGeodesicNode(const SingularityGeodesicNode &vgn)
            : BaseGeodesicNode(vgn),
              vertexID(vgn.vertexID),
              mainDirection(vgn.mainDirection),
              oppositeSinkNode(vgn.oppositeSinkNode),
              exitFace(vgn._exitFace)
        { }

        /**
         * @brief operator = Copy operator.
         * @param vgn
         * @return
         */
        SingularityGeodesicNode & operator =(const SingularityGeodesicNode &vgn) {
            if (this != &vgn) {
                BaseGeodesicNode::operator =(vgn);
                vertexID = vgn.vertexID;
                mainDirection = vgn.mainDirection;
                exitFace = vgn.exitFace;
            }
            return *this;
        }

#if __cplusplus >= 201103L
        /**
         * @brief SingularityGeodesicNode Move constructor.
         * @param vgn
         */
        SingularityGeodesicNode(SingularityGeodesicNode &&vgn)
            : BaseGeodesicNode(std::move(vgn)),
              vertexID(vgn.vertexID),
              mainDirection(vgn.mainDirection),
              oppositeSinkNode(vgn.oppositeSinkNode),
              exitFace(vgn.vgn.exitFace)
        { }

        /**
         * @brief operator = Move operator.
         * @param vgn
         * @return
         */
        SingularityGeodesicNode & operator =(SingularityGeodesicNode &&vgn) {
            if (this != &vgn) {
                BaseGeodesicNode::operator =(std::move(vgn));
                vertexID = vgn.vertexID;
                mainDirection = vgn.mainDirection;
                exitFace=vgn.exitFace;
            }
            return *this;
        }
#endif
    };

#endif //SINGULARITY_NODE_ANISOTROPIC_GEODESIC
