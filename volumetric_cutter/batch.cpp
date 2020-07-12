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

#include "batch.h"
#include "cut.h"
#include "polyhedral_decomposition.h"
#include "subdivision_helper.h"
#include "mesh_extractor.h"
#include "finalization.h"
#include "mesh_smoother.h"
#include <cinolib/export_hexahedra.h>
#include <cinolib/string_utilities.h>

extern std::string output_name;
extern std::string model_name;

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void run_batch(GlobalState & state)
{
    state.profiler.push("LoopyCuts");

    uint cut_count = 0;
    bool converged = false;

    for(uint i=0; i<state.loops.size(); ++i)
    {
        if(state.loops.at(i).type==CONCAVE || !converged)
        {
            if(state.loops.at(i).type==CONVEX) continue;

            cut(state, i);
            state.profiler.push("Mesh Extractor");
            MeshExtractor me(state.m_vol);
            state.m_poly = me.mm;
            converged = me.converged(state.m_vol);
            state.profiler.pop("\n");
            ++cut_count;
        }
    }

    auto model = get_file_name(model_name,false);

    std::cout << cut_count << " cuts performed" << std::endl;
    finalize_block_decomposition(state);
    state.m_poly.save((output_name + "/" + model + "_mm.hedra").c_str());
    SubdivisionHelper sh(state.m_vol, state.m_poly);
    state.m_poly = sh.subdivide();
    state.m_poly.save((output_name + "/" + model + "_mm_subdivided.hedra").c_str());
    smoother(state.m_poly, state.m_srf);
    state.m_poly.save((output_name + "/" + model + "_mm_subdivided_smoothed.hedra").c_str());
    classify_polyhedra(state.m_poly);
    state.profiler.pop();

    // MAKE HEXMESH
    Hexmesh<MM,MV,ME,MF,MP> hm;
    export_hexahedra(state.m_poly, hm);
    uint  nh = hm.num_polys();
    uint  np = state.m_poly.num_polys();
    if(nh==np)
    {
        hm.poly_fix_orientation();
        hm.save((output_name + "/" + model + "_hex.mesh").c_str());
    }
}
