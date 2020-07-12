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
