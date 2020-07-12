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

#include "mesh_smoother.h"

#include <cinolib/smoother.h>
#include <cinolib/export_surface.h>

void smoother(MetaMesh & m, const SrfMesh & target)
{
    std::cout << "smooth meta mesh" << std::endl;
    Polygonmesh<MM,MV,ME,MF> srf;
    std::unordered_map<uint,uint> m2srf_vmap, srf2m_vmap;
    export_surface(m, srf, m2srf_vmap, srf2m_vmap);
    for(uint eid=0; eid<m.num_edges(); ++eid)
    {
        // creases
        if(m.edge_data(eid).flags[MARKED])
        {
            if(!CONTAINS(m2srf_vmap, m.edge_vert_id(eid,0))) continue;
            if(!CONTAINS(m2srf_vmap, m.edge_vert_id(eid,1))) continue;
            int id = srf.edge_id(m2srf_vmap.at(m.edge_vert_id(eid,0)),
                                 m2srf_vmap.at(m.edge_vert_id(eid,1)));
            assert(id>=0);
            srf.edge_data(id).flags[MARKED] = true;
        }
    }
    SmootherOptions opt;
    opt.w_feature = 10000;
    opt.w_corner = 10000;
    opt.n_iters = 100;
    opt.w_laplace = 0.1;
    mesh_smoother(srf, target, opt);
    for(uint vid=0; vid<srf.num_verts(); ++vid)
    {
        m.vert(srf2m_vmap.at(vid)) = srf.vert(vid);
    }

    for(uint i=0; i<10; ++i)
    {
        for(uint vid=0; vid<m.num_verts(); ++vid)
        {
            if(!m.vert_is_on_srf(vid))
            {
                vec3d pos = m.vert(vid);
                for(uint nbr : m.adj_v2v(vid)) pos += m.vert(nbr);
                m.vert(vid) = pos/static_cast<double>(m.adj_v2v(vid).size()+1);
            }
        }
    }
    m.update_normals();
    m.updateGL();
}
