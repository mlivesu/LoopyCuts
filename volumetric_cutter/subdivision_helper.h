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

#ifndef SUBDIVISION_HELPER_H
#define SUBDIVISION_HELPER_H

#include "definitions.h"
#include <unordered_map>
#include <cinolib/octree.h>


//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

typedef struct
{
    Trimesh<> m;                              // mesh of the local patch
    std::vector<uint>  corners;               // patch corners (local vids)
    std::unordered_map<uint,uint> mm2corners; // connects MM verts with local corners
    bool flawed = false;                      // if anything unexpected occurs... discard it
    Octree octree = Octree(4,100);
}
Patch;

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

class SubdivisionHelper
{
    public:

        SubdivisionHelper(const TetMesh  & m,
                          const MetaMesh & mm);

        //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        MetaMesh subdivide();

        //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        void extract_surface_patches();
        void make_uv_maps();
        bool map_to_polygon(Trimesh<> & m, const std::vector<uint> & corners);

        //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    private:

        MetaMesh mm_refined;
        const TetMesh  & m;
        const MetaMesh & mm;

        std::unordered_map<int,Patch> patches; // one for each MM face
};

#endif // SUBDIVISION_HELPER_H
