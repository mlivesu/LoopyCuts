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

#ifndef EDGE_NODE_ANISOTROPIC_GEODESIC
#define EDGE_NODE_ANISOTROPIC_GEODESIC

#include "base_geodesic_node.h"

    /**
     * @brief The EdgeGeodesicNode class
     */
    template <class CoordType>
    class EdgeGeodesicNode : public BaseGeodesicNode<CoordType>
    {
    public:
        typedef BaseGeodesicNode<CoordType> BaseGeodesicNode;

        size_t  faceID;         ///< The index of the face the node belongs to.
        int     edgeIdx;        ///< The index of the edge the node belongs to.
        int     M4Dir;          ///< The number of the vector of the cross field of M4.
        int     OppositeID;     ///< The index of the opposite Node.
        int     OppositeSameFID;///< The index of the opposite Node  considering the same Face
        bool    IsBorder;       ///< true if the node is on border.

        /**
         * @brief EdgeGeodesicNode Constructor.
         * @param _Face
         * @param _edgeIdx
         * @param _M4Dir
         * @param _pos
         */
        EdgeGeodesicNode(const size_t _faceID,
                         const int _edgeIdx,
                         const int _M4Dir,
                         const CoordType &_pos,
                         bool _IsBorder)
            : BaseGeodesicNode(EdgeNodeType, _pos),
              faceID(_faceID),
              edgeIdx(_edgeIdx),
              M4Dir(_M4Dir),
              IsBorder(_IsBorder)
        {
            OppositeID=-1;
            OppositeSameFID=-1;
        }

        /**
         * @brief EdgeGeodesicNode Copy constructor.
         * @param egn
         */
        EdgeGeodesicNode(const EdgeGeodesicNode &egn)
            : BaseGeodesicNode(egn),
              faceID(egn.faceID),
              edgeIdx(egn.edgeIdx),
              M4Dir(egn.M4Dir),
              OppositeID(egn.OppositeID),
              OppositeSameFID(egn.OppositeSameFID),
              IsBorder(egn.IsBorder)
        { }

        /**
         * @brief operator = Copy operator.
         * @param egn
         * @return
         */
        EdgeGeodesicNode & operator =(const EdgeGeodesicNode &egn) {
            if (this != &egn) {
                BaseGeodesicNode::operator =(egn);
                faceID = egn.faceID;
                edgeIdx = egn.edgeIdx;
                M4Dir = egn.M4Dir;
                OppositeID=egn.OppositeID;
                OppositeSameFID=egn.OppositeSameFID;
                IsBorder=egn.IsBorder;
            }
            return *this;
        }

#if __cplusplus >= 201103L
        /**
         * @brief EdgeGeodesicNode Move constructor.
         * @param egn
         */
        EdgeGeodesicNode(EdgeGeodesicNode &&egn)
            : BaseGeodesicNode(std::move(egn)),
              faceID(egn.faceID),
              edgeIdx(egn.edgeIdx),
              M4Dir(egn.M4Dir),
              OppositeID(egn.OppositeID),
              OppositeSameFID(egn.OppositeSameFID),
              IsBorder(egn.IsBorder)
        { }

        /**
         * @brief operator = Move operator.
         * @param egn
         * @return
         */
        EdgeGeodesicNode & operator =(EdgeGeodesicNode &&egn) {
            if (this != &egn) {
                BaseGeodesicNode::operator =(std::move(egn));
                faceID = egn.faceID;
                edgeIdx = egn.edgeIdx;
                M4Dir = egn.M4Dir;
                OppositeID=egn.OppositeID;
                OppositeSameFID=egn.OppositeSameFID;
                IsBorder=egn.IsBorder;
            }
            return *this;
        }
#endif
    };

#endif //EDGE_NODE_ANISOTROPIC_GEODESIC
