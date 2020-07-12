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

#ifndef NEIGH_INFO_ANISOTROPIC_GEODESIC
#define NEIGH_INFO_ANISOTROPIC_GEODESIC

#include <list>
#include <deque>
#include "neigh_info.h"

    /**
     * @brief The NeighborInfo struct
     */
    template <class ScalarType>
    struct NeighborInfo {
        size_t      nodeID;
        ScalarType  angle;
        ScalarType  length;
        bool        Valid;
        ScalarType  Weight;

        /**
         * @brief NeighborInfo Constructor.
         * @param _nodeID
         * @param _angle
         * @param _length
         */
        NeighborInfo(size_t _nodeID,
                     ScalarType _angle,
                     ScalarType _length)
            : nodeID(_nodeID),
              angle(_angle),
              length(_length) {Valid=true;Weight=1;}

        /**
         * @brief NeighborInfo Copy constructor.
         * @param ni
         */
        NeighborInfo(const NeighborInfo &ni)
            : nodeID(ni.nodeID),
              angle(ni.angle),
              length(ni.length),
              Valid(ni.Valid),
              Weight(ni.Weight){ }

        /**
         * @brief operator = Copy operator.
         * @param ni
         * @return
         */
        NeighborInfo & operator =(const NeighborInfo &ni) {
            if (this != &ni) {
                nodeID = ni.nodeID;
                angle = ni.angle;
                length = ni.length;
                Valid = ni.Valid;
                Weight=ni.Weight;
            }
            return *this;
        }

        inline bool operator <(const NeighborInfo &ni) const {
            if (angle != ni.angle)
                return angle < ni.angle;
            if (length != ni.length)
                return length < ni.length;
            return nodeID < ni.nodeID;
        }

        inline bool operator ==(const NeighborInfo &ni) const {
            return nodeID == ni.nodeID && angle == ni.angle && length == ni.length && Valid==ni.Valid;
        }

        /**
         * @brief geoDistance computes the geodesic distance to this neighbor.
         * @param max_angle
         * @param drift_penalty
         * @return
         */
        inline ScalarType geoDistance(const ScalarType &max_angle = 45,
                                      const ScalarType &drift_penalty = 2)
        {
            if (max_angle == ScalarType(0))
                return (length*Weight);
            //return length * (ScalarType(1) + drift_penalty * angle / max_angle)*Weight;
            ScalarType drift=pow(angle/max_angle,2);
            //ScalarType drift=angle/max_angle;
//            if (Weight<1)
//                std::cout<<"W "<<Weight<<std::endl;
//            else
//                std::cout<<"W 1"<<Weight<<std::endl;
            return length * (ScalarType(1) + drift_penalty * drift)*Weight;
        }
    };

//    typedef std::list<NeighborInfo>                 NeighborsList;
//    typedef typename NeighborsList::iterator        NeighborsIterator;
//    typedef typename NeighborsList::const_iterator  NeighborsConstIterator;

#endif //NEIGH_INFO_ANISOTROPIC_GEODESIC
