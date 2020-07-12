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

#include "finalization.h"
#include "mesh_extractor.h"
#include "cut.h"

extern bool batch_mode;

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

std::vector<uint> edge_loop_around_vert(const MetaMesh                 & m,
                                        const uint                       pid,
                                        const uint                       vid,
                                        const std::unordered_set<uint> & new_verts)
{
    assert(m.poly_contains_vert(pid,vid));

    auto seed = new_verts.begin();
    while(!m.poly_contains_vert(pid,*seed)) ++seed;
    assert(m.poly_contains_vert(pid,*seed));

    // trace edge loop around pid
    std::vector<uint> loop = { *seed };
    std::unordered_set<uint> visited;
    visited.insert(*seed);
    while(true)
    {
        bool found_next = false;
        for(uint id : m.poly_v2v(pid, loop.back()))
        {
            if(CONTAINS(new_verts,id) && DOES_NOT_CONTAIN(visited,id))
            {
                assert(m.edge_id(loop.back(),id)>=0);
                visited.insert(id);
                loop.push_back(id);
                found_next = true;
                break;
            }
        }
        if(!found_next) break;
        assert(loop.size()<=new_verts.size());
    }
    return loop;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void remove_srf_val4_vertex(MetaMesh & m, const uint vid)
{
    if(!m.vert_is_on_srf(vid))   return;
    if(m.adj_v2p(vid).size()>1)  return;
    if(m.adj_v2v(vid).size()!=4) return;
    assert(m.poly_vert_valence(m.adj_v2p(vid).front(),vid)==4);

    std::cout << "Remove val 4 vert " << vid << " from meta mesh" << std::endl;

    uint nv = m.num_verts();

    // split edges
    std::vector<ipair> e_to_split;
    for(uint eid : m.adj_v2e(vid))
    {
        e_to_split.push_back(std::make_pair(m.edge_vert_id(eid,0),
                                            m.edge_vert_id(eid,1)));
    }
    std::vector<uint> crease_points;
    for(auto e : e_to_split)
    {
        int  eid    = m.edge_id(e);
        bool crease = m.edge_data(eid).flags[MARKED];
        assert(eid>=0);
        double t = (m.edge_vert_id(eid,0)==vid)?0.1:0.9;
        uint new_vid = m.edge_split(eid, m.edge_sample_at(eid, t));
        if(crease)
        {
            for(uint eid : m.adj_v2e(new_vid)) m.edge_data(eid).flags[MARKED] = true;
            crease_points.push_back(new_vid);
        }
    }

    // split faces
    std::vector<std::vector<uint>> f_to_split;
    for(uint fid : m.adj_v2f(vid))
    {
        f_to_split.push_back(m.face_verts_id(fid));
    }
    std::unordered_set<uint> new_verts;
    for(auto f : f_to_split)
    {
        int fid = m.face_id(f);
        assert(fid>=0);
        assert(m.face_contains_vert(fid,vid));

        std::vector<uint> endpoints;
        for(uint vid : f) if(vid>=nv) endpoints.push_back(vid);
        assert(endpoints.size()==2);

        new_verts.insert(endpoints.front());
        new_verts.insert(endpoints.back());
        m.face_split_along_new_edge(fid, endpoints.front(), endpoints.back());
    }
    assert(new_verts.size()==e_to_split.size());

    // update polys
    std::vector<uint> loop = edge_loop_around_vert(m, m.adj_v2p(vid).front(), vid, new_verts);
    uint new_face = m.face_add(loop);

    if(crease_points.size()==2) // mark sharp features again
    {
        for(uint vid : new_verts)
        {
            if(m.verts_are_adjacent(vid, crease_points.front()) &&
               m.verts_are_adjacent(vid, crease_points.back()))
            {
                int e0 = m.edge_id(vid, crease_points.front());
                int e1 = m.edge_id(vid, crease_points.back());
                assert(e0 >= 0);
                assert(e1 >= 0);
                m.edge_data(e0).flags[MARKED] = true;
                m.edge_data(e1).flags[MARKED] = true;
                break;
            }
        }
    }

    std::vector<std::vector<uint>> p_to_update;
    for(uint pid : m.adj_v2p(vid))
    {
        p_to_update.push_back(m.poly_faces_id(pid));
    }
    for(auto p : p_to_update)
    {
        int pid = m.poly_id(p);
        assert(pid>=0);
        std::vector<uint> f;
        std::vector<bool> w;
        for(uint fid : m.poly_faces_id(pid))
        {
            if(!m.face_contains_vert(fid,vid))
            {
                f.push_back(fid);
                w.push_back(m.poly_face_winding(pid,fid));
            }
        }
        f.push_back(new_face);
        w.push_back(true); // TODO: fix winding!
        m.poly_add(f,w);
    }
    m.vert_remove(vid);
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

// - adds cuts that fix valence 2 verts
// - remove cuts that generated non fixable valence 2 verts
// - removes valence 4 verts that may be generated by the previous two steps
void finalize_block_decomposition(GlobalState & state)
{
    auto M     = state.m_vol;
    auto MM    = state.m_poly;
    bool fail  = false;
    uint count = 0;

    while(true)
    {
        MeshExtractor me(state.m_vol);
        if(++count>5 || me.mm.num_polys()==0) // *|| !me.converged(state.m_vol)*/)
        {
            fail = true;
            break;
        }

        std::unordered_set<uint> cuts_to_do;    // cuts that would fix some existing val2 verts
        std::unordered_set<uint> cuts_to_undo;  // cuts that created val2 verts that cannot be removed
        me.post_convergence_cut_analysis(state.m_vol, state.loops, cuts_to_do, cuts_to_undo);

        //// this addresses valence two verts that arise when a concave feature is crossed by
        //// a loop that defines a surface edge on the surface. Basically, the two cuts along the
        //// concavity define three elements, one of which has only edges on the surface, and not faces.
        //// such element will have the intersection vertex, with valence two
        //// (ORRIBLE) current solution: don't cut twice along the concavity
        //for(uint vid=0; vid<me.mm.num_verts(); ++vid)
        //{
        //    bool three_incident_polys = me.mm.adj_v2p(vid).size()==3;
        //    bool val2_wrt_some_poly   = false;
        //    for(uint pid : me.mm.adj_v2p(vid))
        //    {
        //        if(me.mm.poly_vert_valence(pid,vid)==2)
        //        {
        //            val2_wrt_some_poly = true;
        //            break;
        //        }
        //    }

        //    if(three_incident_polys && val2_wrt_some_poly)
        //    {
        //        std::set<uint> cuts;
        //        uint m_vid = me.mm.vert_data(vid).vid_on_m;
        //        for(uint fid : state.m_vol.adj_v2f(m_vid))
        //        {
        //            if(state.m_vol.face_data(fid).cut_id>=0) cuts.insert(state.m_vol.face_data(fid).cut_id);
        //        }
        //        if(cuts.size()==2)
        //        {
        //            auto c0 = cuts.begin();
        //            auto c1 = cuts.rbegin();
        //            if(DOES_NOT_CONTAIN(cuts_to_undo, *c0) && DOES_NOT_CONTAIN(cuts_to_undo, *c1))
        //            {
        //                cuts_to_undo.insert(*c0);
        //            }
        //        }
        //    }
        //}

        if(cuts_to_do.empty() && cuts_to_undo.empty()) break; // converged

        PRINT(cuts_to_do,   "CUTS TO DO");
        PRINT(cuts_to_undo, "CUTS TO UNDO");

        for(uint lid : cuts_to_undo)
        {
            revert_cut(state, lid);
        }

        for(uint lid : cuts_to_do)
        {
            Loop    & l = state.loops.at(lid);
            TetMesh & m = state.m_vol;

            // Jan 17 2020: try this heuristic: if the mesh is smooth, add the loop as MM edges,
            // if it contains concave creases, then try cutting
            // this is because loops as edges make val 2 verts when intersect concave crease cut through twice
            if(state.loops.front().type==CONCAVE) cut(state, lid);
            else l.reverted = true;

            // if cutting along te loop was impossible, make a crease out of it
            if(l.reverted)
            {
                l.type = TOP_RELEVANT;
                std::vector<uint>  vids;        // IDS of loop points
                std::vector<vec3d> points;      // coordinates of loop points
                std::vector<vec3d> srf_normals; // surface normals at each loop point
                std::vector<vec3d> cut_normals; // cut normal at each loop point
                state.loops.trace_loop(lid, vids, points, srf_normals, cut_normals);

                for(auto i=vids.begin(),j=i+1; j<vids.end(); ++i,++j)
                {
                    int eid = m.edge_id(*i,*j);
                    assert(eid>=0);
                    m.edge_data(eid).topologically_needed = true; // makes sure these edges will form edges in the meta mesh
                }
                if(l.closed)
                {
                    int eid = m.edge_id(*vids.begin(),*vids.rbegin());
                    assert(eid>=0);
                    m.edge_data(eid).topologically_needed = true;
                }
            }
        }
    }

    // If I failed to make a full hexmesh, revert
    MeshExtractor me(state.m_vol);
    for(uint pid=0; pid<me.mm.num_polys(); ++pid)
    {
        if(!me.mm.poly_is_hexable_w_midpoint(pid))
        {
            fail = true;
            break;
        }
    }

    if(!fail)
    {
        state.m_poly = me.mm;
        if(!me.converged(state.m_vol)) std::cout << "\n\n\nWARNING: Non hexa has been removed, but geometric convergence is now lost...by how much?\n\n\n" << std::endl;
    }
    else
    {
        std::cout << "\n\n\nWARNING: Any attempt to make a full hexa mesh failed. Revert to converged hybrid mesh...\n\n\n" << std::endl;
        state.m_poly = MM;
        state.m_vol  = M;
    }

    // loops transformed into creases that intersect will generte valence 4 verts.
    // here I remove them all
    for(uint vid=0; vid<state.m_poly.num_verts(); ++vid)
    {
        for(uint pid : state.m_poly.adj_v2p(vid))
        {
            if(state.m_poly.poly_vert_valence(pid,vid)==4) remove_srf_val4_vertex(state.m_poly, vid);
        }
    }

    if(!batch_mode) state.m_poly.updateGL();
}
