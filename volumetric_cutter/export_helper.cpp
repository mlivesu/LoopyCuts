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

#include "export_helper.h"
#include <cinolib/export_hexahedra.h>

void exporter(const GlobalState & state, const std::string filename)
{
    std::unordered_map<uint,uint> v_map;
    Hexmesh<MM,MV,ME,MF,MP> hm;
    export_hexahedra(state.m_poly, hm, v_map);
    if(hm.num_polys()>0)
    {
        hm.poly_fix_orientation(); // meta mesh does not have globally consistent winding numbers...
        hm.save((filename + "_hex.mesh").c_str());
    }
    state.m_vol.save((filename + "_tet.mesh").c_str());

    if(state.m_poly.num_polys() > hm.num_polys())
    {
        state.m_poly.save((filename + "_midpoint.hedra").c_str());

        FILE *f = fopen((filename + "_vmap.txt").c_str(), "w");
        assert(f);
        for(auto obj : v_map) fprintf(f,"%d %d\n", obj.first, obj.second);
        fclose(f);
    }
}
