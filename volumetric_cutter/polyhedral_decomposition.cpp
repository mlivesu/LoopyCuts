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

#include <unordered_map>
#include <cinolib/ANSI_color_codes.h>
#include "mesh_extractor.h"
#include "polyhedral_decomposition.h"

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void polyhedron_details(const MetaMesh & m, const uint pid)
{
    std::cout << "Polyhedron #" << pid << ": " << m.verts_per_poly(pid) << "v/" << m.edges_per_poly(pid) << "e/" << m.faces_per_poly(pid) << "f" << std::endl;
    std::unordered_map<uint,uint> valence_map;
    for(uint vid : m.adj_p2v(pid))
    {
        ++valence_map[m.poly_vert_valence(pid,vid)];
    }
    for(auto obj : valence_map)
    {
        std::cout << "\t" << obj.second << " valence " << obj.first << std::endl;
    }
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void classify_polyhedra(MetaMesh & m)
{
    uint n_hexa    = 0;
    uint n_prism   = 0;
    uint n_hexable = 0;
    uint n_other   = 0;
    for(uint pid=0; pid<m.num_polys(); ++pid)
    {
        if(m.poly_is_hexahedron(pid))
        {
            ++n_hexa;
            m.poly_data(pid).label = 0;
            m.poly_data(pid).color = Color::PASTEL_GREEN();
        }
        else if(m.poly_is_prism(pid))
        {
            ++n_prism;
            m.poly_data(pid).label = 1;
            m.poly_data(pid).color = Color::PASTEL_VIOLET();
        }
        else if(m.poly_is_hexable_w_midpoint(pid))
        {
            ++n_hexable;
            m.poly_data(pid).label = 2;
            m.poly_data(pid).color = Color::PASTEL_ORANGE();
        }
        else
        {
            ++n_other;
            m.poly_data(pid).label = 3;
            m.poly_data(pid).color = Color::PASTEL_RED();
        }
    }
    std::cout << ANSI_fg_color_green   << n_hexa    << " hexa - "
              << ANSI_fg_color_blue    << n_prism   << " prisms - "
              << ANSI_fg_color_red     << n_hexable << " hexable with midpoint - "
              << ANSI_fg_color_red     << n_other   << " others"
              << ANSI_fg_color_default << std::endl;
    m.updateGL();
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void extract_polyhedral_mesh(GlobalState & state)
{
    state.profiler.push("update meta mesh");
    MeshExtractor mesh_ex(state.m_vol);
    mesh_ex.converged(state.m_vol);
    state.profiler.pop();
    state.m_poly = mesh_ex.mm;
    classify_polyhedra(state.m_poly);
    std::cout << std::endl;
}
