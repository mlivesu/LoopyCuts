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

#ifndef DEFINITIONS_H
#define DEFINITIONS_H

#include <cinolib/meshes/meshes.h>

using namespace cinolib;

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

// loop types
enum
{
    CONCAVE,      // loop on concave crease
    CONVEX,       // loop on convex crease
    REGULAR,      // loop not on a crease
    TOP_RELEVANT, // loop that must appear as surface edge in the meta mesh (to balance vert valences)
    NONE,         // not a loop
};

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

static std::vector<std::string> TYPES_STR =
{
    "CONCAVE",
    "CONVEX",
    "REGULAR",
    "NONE",
};

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

typedef struct : Edge_std_attributes
{
    bool  sharp_crease         = false;
    bool  topologically_needed = false;
    int   loop_ids[2]          = { -1, -1 };
    vec3d srf_normals[2];
    bool  has_two_loops() const { return loop_ids[0]>=0 && loop_ids[1]>=0; }
    bool  must_be_mm_edge() const { return sharp_crease || topologically_needed; }
    bool  belongs_to_loop(const uint id) const { return loop_ids[0]==(int)id || (has_two_loops() && loop_ids[1]==(int)id); }
    vec3d srf_normal(const uint id) const { assert(belongs_to_loop(id)); if(loop_ids[0]==(int)id) return srf_normals[0]; else return srf_normals[1]; }
    int   type = NONE;  // this I don't know
}
E;

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

// SURFACE and VOLUME MESH attributes
typedef          Mesh_std_attributes                                              M;
typedef struct : Vert_std_attributes       { bool on_loop = false; }              V;
typedef struct : Polygon_std_attributes    { int cut_id = -1; }                   F;
typedef          Polyhedron_std_attributes                                        P;

typedef DrawableTrimesh<M,V,E,F>   SrfMesh;
typedef DrawableTetmesh<M,V,E,F,P> TetMesh;

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

// META MESH attributes
typedef          Mesh_std_attributes                                                                                                  MM;
typedef struct : Vert_std_attributes       { uint vid_on_m;                                                                         } MV;
typedef struct : Edge_std_attributes       { std::vector<uint> v_chain_on_m; std::unordered_set<uint> face_labels; int type = NONE; } ME;
typedef struct : Polygon_std_attributes    {                                                                                        } MF;
typedef struct : Polyhedron_std_attributes {                                                                                        } MP;

typedef DrawablePolyhedralmesh<MM,MV,ME,MF,MP> MetaMesh;

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

typedef struct
{
    std::vector<vec3d>       drawlist;           // serialized list of edges (for rendering)
    uint                     starting_vert;      // id of the vert from which start tracing the loop (any if it's a closed loop)
    std::vector<uint>        verts;              // ordered vert list (ids of both m_srf && m_vol)
    std::unordered_set<uint> polys;              // used for manual picking
    int                      type;               // either CONCAVE, CONVEX or REGULAR
    bool                     closed;             // true if it is a closed loop
    bool                     flawed;             // the loop has two consecutive intersections with some other loop. Cutting along it will produce non hexa elements
    bool                     used = false;       // yes if it has already been used in a volumetric cut
    bool                     reverted = false;   // yes if the loop created e.g. some bubble, hence was computed but not used
    bool                     srf_bubble = false; // true if the loop traverses the surface, creating a new loop
    bool                     active = false;     // used for visualization
    bool                     Nico_bug = false;   // some convex loops seem to be duplicated, I am skipping them...
}
Loop;

#endif // DEFINITIONS_H
