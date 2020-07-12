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

#ifndef MESH_EXTRACTOR_H
#define MESH_EXTRACTOR_H

#include "definitions.h"
#include "loops.h"
#include "state.h"
#include <unordered_map>

class MeshExtractor
{
    public:

        MeshExtractor(TetMesh & m);

        //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        MetaMesh mm;

        //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        void init      (TetMesh & m);
        void make_verts(TetMesh & m);
        void make_edges(TetMesh & m); // new verts can be added (and marked in m)
        void make_faces();
        void make_polys(TetMesh & m);

        //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        bool converged              (const TetMesh & m) const;
        bool geometric_convergence  (const TetMesh & m) const; // TODO: check edge/arc length mismatch
        bool topological_convergence(const TetMesh & m) const;

        //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        void post_convergence_cut_analysis(const TetMesh            & m,
                                           const Loops              & loops,
                                           std::unordered_set<uint> & cuts_to_do,     // cuts that would fix some existing val2 verts
                                           std::unordered_set<uint> & cuts_to_undo);  // cuts that created val2 verts that cannot be removed

        //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    protected:

        std::unordered_map<uint,uint> m2mm_vmap; // maps vert IDs of the tetmesh into vert IDs in the meta mesh
        std::unordered_map<uint,uint> m2mm_fmap; // maps per face labels of the tetmesh into face IDs in the meta mesh

        //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        uint polys_non_manifold    = 0;
        uint polys_with_high_genus = 0;

        //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        std::map<ipair,std::vector<std::vector<uint>>> buggy_chains;

        //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        std::vector<uint>        make_edge(const TetMesh & m, const uint start_v, const uint start_e) const;
        std::vector<uint>        make_face(const uint start_e, const int label) const;
        std::unordered_set<uint> make_poly(TetMesh &m, const uint start_p) const;
        bool                     remove_degenerate_edges(); // kill vertices with valence 1 and 2
        void                     update_v_map();
        void                     update_f_map();
        std::unordered_set<uint> face_labels_incident_to_edge(const TetMesh & m, const uint mm_eid);
        bool                     topological_checks(const std::unordered_set<uint> & poly);

        //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        std::unordered_set<uint> srf_face_labels;
};

#endif // MESH_EXTRACTOR_H
