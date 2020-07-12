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

#include "subdivision_helper.h"
#include <cinolib/geometry/n_sided_poygon.h>
#include <cinolib/sampling.h>
#include <cinolib/harmonic_map.h>
#include <cinolib/subdivision_midpoint.h>

SubdivisionHelper::SubdivisionHelper(const TetMesh & m, const MetaMesh & mm)
: m(m)
, mm(mm)
{
    extract_surface_patches();
    make_uv_maps();
    subdivide();
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void SubdivisionHelper::extract_surface_patches()
{
    /* meta mesh extrctor labels all facets in the tetmesh and associates
     * a unique label to each meta mesh face, creating a connection between
     * the two face sets. Here I exploit this connection to retrieve the patches
     * to be uv-maped
    */

    for(uint mm_fid=0; mm_fid<mm.num_faces(); ++mm_fid)
    {
        if(!mm.face_is_on_srf(mm_fid)) continue;
        int l = mm.face_data(mm_fid).label;
        if(l<0) continue; // from the second run on, no uv maps available (will do standard midpoint)

        std::cout << "extracting patch for MM surface face " << mm_fid << ". Face label on M is " << l << std::endl;

        Patch p;
        std::unordered_map<uint,uint> vmap;
        for(uint m_fid=0; m_fid<m.num_faces(); ++m_fid)
        {
            if(m.face_data(m_fid).label==l)
            {
                uint tri[3];
                for(int i=0; i<3; ++i)
                {
                    uint vid = m.face_vert_id(m_fid,i);
                    auto query = vmap.find(vid);
                    if(query==vmap.end())
                    {
                        uint fresh_id = p.m.vert_add(m.vert(vid));
                        vmap[vid] = fresh_id;
                        tri[i] = fresh_id;
                    }
                    else tri[i] = query->second;
                }
                p.m.poly_add(tri[0], tri[1], tri[2]);
            }
        }

        for(uint mm_corner : mm.adj_f2v(mm_fid))
        {
            uint m_corner = mm.vert_data(mm_corner).vid_on_m;
            uint p_corner = vmap.at(m_corner);
            p.corners.push_back(p_corner);
            p.mm2corners[mm_corner] = p_corner;
        }
        patches[l] = p;
        //patches[l].m.save(("/Users/cino/Desktop/" + std::to_string(l) + ".obj").c_str());
    }
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void SubdivisionHelper::make_uv_maps()
{
    for(auto & p : patches)
    {
        std::cout << "uv mapping patch " << p.first << "(corners " << p.second.corners.size() << ")" << std::endl;
        if(map_to_polygon(p.second.m, p.second.corners)) p.second.octree.build_from_mesh_polys(p.second.m);
        else
        {
            p.second.flawed = true;
            p.second.m.save(("/Users/cino/Desktop/errors/" + std::to_string(p.first) + ".obj").c_str());
        }
        //for(auto c : p.second.corners_uv) std::cout << "\t" << c.first << " => " << c.second << std::endl;
    }
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

bool SubdivisionHelper::map_to_polygon(Trimesh<> & m, const std::vector<uint> & corners)
{
    // make sure it's a topological disk
    //if(m.Euler_characteristic()!=1)
    //{
    //    std::cout << "WARNING: patch is not a topological disk!" << std::endl;
    //    return false;
    //}

    // get the ordered list of boundary vertices, starting from a corner
    m.vert_unmark_all();
    for(uint vid : corners) m.vert_data(vid).flags[MARKED] = true;
    std::vector<uint> border = m.get_ordered_boundary_vertices();
    if(border.empty())
    {
        std::cout << "WARNING: non orientable patch boundary!" << std::endl;
        return false;
    }
    CIRCULAR_SHIFT_VEC(border, corners.front());

    // split the boundary into n edges, with n = #corners
    std::vector<std::vector<uint>> edges;
    for(uint i=0; i<border.size(); ++i)
    {
        std::vector<uint> e = { border.at(i) };
        for(uint j=i+1; j<border.size() && !m.vert_data(border.at(j)).flags[MARKED]; ++j,++i)
        {
            e.push_back(border.at(j));
        }
        e.push_back(border.at((i+1)%border.size()));
        edges.push_back(e);
    }

    std::vector<vec3d> poly = n_sided_polygon(vec3d(0,0,0), corners.size(), 1.0);
    std::map<uint,vec3d> dirichlet_bcs;
    for(uint i=0; i<poly.size(); ++i)
    {
        std::vector<vec3d> e_bcs = sample_within_interval(poly.at(i), poly.at((i+1)%poly.size()), edges.at(i).size());
        for(uint j=0; j<e_bcs.size(); ++j) dirichlet_bcs[edges.at(i).at(j)] = e_bcs.at(j);
    }

    // map the interior vertices
    m.copy_xyz_to_uvw(UVW_param); // save XYZ coordinates in UVW
    m.vector_verts() = harmonic_map_3d(m, dirichlet_bcs, 1, COTANGENT); // map to disk

    return true;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

MetaMesh SubdivisionHelper::subdivide()
{
    // do standard midpoint subdivision
    std::unordered_map<uint,uint> edge_verts;
    std::unordered_map<uint,uint> face_verts;
    std::unordered_map<uint,uint> poly_verts;
    subdivision_midpoint(mm, mm_refined, edge_verts, face_verts, poly_verts);

    // restore sharp creases;
    for(auto e : edge_verts)
    {
        uint mm_eid = e.first;
        bool crease = mm.edge_data(mm_eid).flags[MARKED];
        if(crease)
        {
            uint vid   = e.second;
            uint mm_v0 = mm.edge_vert_id(mm_eid,0);
            uint mm_v1 = mm.edge_vert_id(mm_eid,1);
            int  e0    = mm_refined.edge_id(mm_v0,vid); assert(e0>=0);
            int  e1    = mm_refined.edge_id(mm_v1,vid); assert(e1>=0);
            mm_refined.edge_data(e0).flags[MARKED]= true;
            mm_refined.edge_data(e1).flags[MARKED]= true;
        }
    }

    // if per patch UV maps are available, reposition surface midpoints using the map
    if(!patches.empty())
    {
        // reposition face midpoints using the uv map of its associated MM face
        for(auto f : face_verts)
        {
            uint mm_fid = f.first;
            if(!mm.face_is_on_srf(mm_fid)) continue;
            uint vid      = f.second;
            int  patch_id = mm.face_data(mm_fid).label;
            Patch & p     = patches.at(patch_id);
            if(p.flawed)
            {
                mm_refined.vert_data(vid).color = Color::RED();
                continue;
            }

            vec3d query(0,0,0);
            for(auto c : p.corners) query += p.m.vert(c);
            query /= static_cast<double>(p.corners.size());

            uint   pid;
            vec3d  pos;
            double dist;
            p.octree.closest_point(query, pid, pos, dist);
            //std::cout << "found face midpoint in uv: " << pos << "\t" << dist << std::endl;
            double bary[3];
            p.m.poly_bary_coords(pid, pos, bary);
            double u = p.m.poly_sample_param_at(pid, bary, U_param);
            double v = p.m.poly_sample_param_at(pid, bary, V_param);
            double w = p.m.poly_sample_param_at(pid, bary, W_param);
            vec3d new_pos(u,v,w);
            mm_refined.vert(vid) = new_pos;
        }

        // reposition edge midpoints using one of the uv map of its incident MM faces
        for(auto e : edge_verts)
        {
            uint mm_eid = e.first;
            if(!mm.edge_is_on_srf(mm_eid)) continue;
            uint vid      = e.second;
            assert(mm.edge_adj_srf_faces(mm_eid).size()>0);
            uint mm_fid   = mm.edge_adj_srf_faces(mm_eid).front();
            int  patch_id = mm.face_data(mm_fid).label;
            Patch & p     = patches.at(patch_id);
            if(p.flawed)
            {
                mm_refined.vert_data(vid).color = Color::RED();
                continue;
            }

            uint mm_v0    = mm.edge_vert_id(mm_eid,0);
            uint mm_v1    = mm.edge_vert_id(mm_eid,1);
            uint c0       = p.mm2corners.at(mm_v0);
            uint c1       = p.mm2corners.at(mm_v1);
            vec3d  query = (p.m.vert(c0) + p.m.vert(c1))*0.5;

            uint   pid;
            vec3d  pos;
            double dist;
            p.octree.closest_point(query, pid, pos, dist);
            //std::cout << "found edge midpoint in uv: " << pos << "\tpid: " << pid << "\tdist: " << dist << std::endl;
            double bary[3];
            p.m.poly_bary_coords(pid, pos, bary);
            double u = p.m.poly_sample_param_at(pid, bary, U_param);
            double v = p.m.poly_sample_param_at(pid, bary, V_param);
            double w = p.m.poly_sample_param_at(pid, bary, W_param);
            vec3d new_pos(u,v,w);
            mm_refined.vert(vid) = new_pos;
        }
    }

    return mm_refined;
}
