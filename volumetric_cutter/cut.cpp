#include "cut.h"
#include <cinolib/interval.h>
#include <cinolib/RBF_Hermite.h>
#include <cinolib/RBF_kernels.h>
#include <cinolib/deg_rad.h>
#include <cinolib/memory_usage.h>
#include "labeling.h"

extern bool batch_mode;

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void tessellate_connected_component_of_isosurface(      TetMesh            & m,
                                                        Loop               & l,
                                                  const std::vector<uint>  & seeds,
                                                  const std::vector<vec3d> & normals,
                                                  const double               iso_value)
{
    // [OVERKILL] split all edges travesed by the isosurface, as if
    // the whole level set should be gloally embedded in the mesh
    typedef std::pair<uint,double> split_data;
    std::set<split_data,std::greater<split_data>> edges_to_split; // from highest to lowest id
    for(uint eid=0; eid<m.num_edges(); ++eid)
    {
        double f0 = m.vert_data(m.edge_vert_id(eid,0)).uvw[0];
        double f1 = m.vert_data(m.edge_vert_id(eid,1)).uvw[0];

        if(is_into_interval<double>(iso_value, f0, f1))
        {
            double alpha = std::fabs(iso_value - f0)/fabs(f1 - f0);
            edges_to_split.insert(std::make_pair(eid,alpha));
            if(m.edge_is_on_srf(eid)) l.srf_bubble = true;
        }
    }
    for(auto e : edges_to_split)
    {
        uint vid = m.edge_split(e.first, e.second);
        m.vert_data(vid).uvw[0] = iso_value;
    }

    // start from the verts in the outer loop (seeds), and propagate along the level set
    // in order to find a connected component of it that realizes the wanted cut.
    // NOTE: in most of the cases each edge will have only one incident face that belongs
    // to the cut. In some unfortunate case there will be two, and only one of them will
    // be on the right side of the cut. The code below correctly selects the right face
    // as the one best aligned with the surface anti normal
    std::queue<uint> q;
    m.vert_unmark_all();
    for(uint vid : seeds) m.vert_data(vid).flags[MARKED] = true;
    for(uint i=0; i<seeds.size(); ++i)
    {
        uint curr = seeds.at(i);
        uint next = seeds.at((i+1)%seeds.size());
         int eid  = m.edge_id(curr, next);
        //assert(eid>=0);

        // with loop pairing there may be multiple loops serialized in the same vectors,
        // thus creating a discontinuity
        if(eid==-1) continue;

        std::vector<uint> v_opp;
        for(uint fid : m.adj_e2f(eid))
        {
            uint vid = m.face_vert_opposite_to(fid, eid);
            if(m.vert_data(vid).uvw[0] == iso_value) v_opp.push_back(vid);
        }

        if(v_opp.size()==1)
        {
            m.vert_data(v_opp.front()).flags[MARKED] = true;
            q.push(v_opp.front());
        }
        else if(v_opp.size()==2)
        {
            vec3d n = normals.at(i);
            vec3d u = m.vert(v_opp.front()) - m.edge_vert(eid,0);
            vec3d v = m.vert(v_opp.back())  - m.edge_vert(eid,0);
            if(u.angle_deg(-n) < v.angle_deg(-n))
            {
                m.vert_data(v_opp.front()).flags[MARKED] = true;
                q.push(v_opp.front());
            }
            else
            {
                m.vert_data(v_opp.back()).flags[MARKED] = true;
                q.push(v_opp.back());
            }
        }
        // not sure why this would happen, but it happens (numerical roundoffs in the cut tracing?)
        else if(v_opp.size()!=1 && v_opp.size()!=2)
        {
            //std::cout << "\n\nSOMETHING WEIRD HAPPENED" << std::endl;
            //for(uint v : v_opp) std::cout << "\t" << v << std::endl;
            continue;
        }
    }
    while(!q.empty())
    {
        uint vid = q.front();
        q.pop();

        assert(m.vert_data(vid).uvw[0] == iso_value);
        for(uint nbr : m.adj_v2v(vid))
        {
            if(!m.vert_data(nbr).flags[MARKED] && m.vert_data(nbr).uvw[0] == iso_value)
            {
                m.vert_data(nbr).flags[MARKED] = true;
                q.push(nbr);
            }
        }
    }
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void find_mates(      GlobalState           & state,
                const Hermite_RBF<CubicRBF> & HRBF,
                      std::vector<uint>     & vids,
                      std::vector<vec3d>    & points,
                      std::vector<vec3d>    & srf_normals,
                      std::vector<vec3d>    & cut_normals)
{
    for(uint lid=0; lid<state.loops.size(); ++lid)
    {
        if(state.loops.at(lid).used || state.loops.at(lid).type==CONVEX) continue;
        double max_dist = -inf_double;
        double max_ang  = -inf_double;

        std::vector<uint>  tmp_vids;
        std::vector<vec3d> tmp_points;
        std::vector<vec3d> tmp_srf_normals;
        std::vector<vec3d> tmp_cut_normals;
        state.loops.trace_loop(lid, tmp_vids, tmp_points, tmp_srf_normals, tmp_cut_normals);

        for(uint i=0; i<tmp_vids.size(); ++i)
        {
            vec3d p    = tmp_points.at(i);
            vec3d n    = tmp_cut_normals.at(i);
            vec3d grad = HRBF.eval_grad(p);
            grad.normalize();
            double dot = grad.dot(n);
            if(dot<0)
            {
                tmp_cut_normals.at(i) = -tmp_cut_normals.at(i);
                dot = -dot;
            }
            double ang = to_deg(acos(dot));
            max_dist = std::max(max_dist, std::fabs(HRBF.eval(p)));
            max_ang  = std::max(max_ang,  ang);
        }
        if(max_dist<0.1 && max_ang<10)
        {
            state.loops.at(lid).used = true;
            if(state.loops.at(lid).type==CONCAVE)
            {
                std::copy(tmp_vids.begin(),        tmp_vids.end(),        std::back_inserter(vids));
                std::copy(tmp_points.begin(),      tmp_points.end(),      std::back_inserter(points));
                std::copy(tmp_srf_normals.begin(), tmp_srf_normals.end(), std::back_inserter(srf_normals));
                std::copy(tmp_cut_normals.begin(), tmp_cut_normals.end(), std::back_inserter(cut_normals));
                std::cout << "loop " << lid << " is interpolated by current cut. Add it to BCs" << std::endl;
            }
            else
            {
                std::cout << "loop " << lid << " is interpolated by current cut. Discard it" << std::endl;
            }
        }
    }
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

bool validate_cut_topology(TetMesh                        & m,
                           const uint                     & cut_id,
                           const std::unordered_set<uint> & faces_on_cut,
                           const std::unordered_set<uint> & edges_on_cut,
                           const std::unordered_set<uint> & /*verts_on_cut*/)
{
    // check faces
    for(uint fid : faces_on_cut)
    {
        if(m.face_data(fid).cut_id>=0)
        {
            std::cout << "Face " << fid << " already belongs to cut " << m.face_data(fid).cut_id << ". Revert!" << std::endl;
            return false;
        }
        if(m.face_is_on_srf(fid))
        {
            std::cout << "Face " << fid << " is on surface: tangent cut. Revert!" << std::endl;
            return false;
        }
        m.face_data(fid).cut_id = cut_id;
    }

    // check edges
    for(uint eid : edges_on_cut)
    {
        uint valence = 0;
        std::unordered_set<uint> e_cuts;
        for(uint fid : m.adj_e2f(eid))
        {
            if(m.face_data(fid).cut_id>=0)
            {
                ++valence;
                e_cuts.insert(m.face_data(fid).cut_id);
            }
        }
        bool vol_regular = (!m.edge_is_on_srf(eid) && valence==2 && e_cuts.size()==1);
        bool mm_srf_edge = ( m.edge_is_on_srf(eid) && valence==1 && e_cuts.size()==1) ||
                           ( m.edge_is_on_srf(eid) && valence==2 && e_cuts.size()==2); // e.g. two cuts along a concave crease
        bool mm_vol_edge = (valence==4 && e_cuts.size()==2);

        if(!vol_regular && !mm_srf_edge && !mm_vol_edge)
        {
            //std::cout << "Edge " << m.edge_vert_ids(eid) << " is neither regular nor a valid mm_edge. (#incident face cuts: " << valence << ", #cuts: " << e_cuts.size() << ", srf: " << m.edge_is_on_srf(eid) << "). Discard!" << std::endl;
            return false;
        }
    }

    // check verts
    //for(uint vid : verts_on_cut)
    //{
    //    uint valence = 0;
    //    std::unordered_set<uint> v_cuts;
    //    for(uint eid : m.adj_v2e(vid))
    //    {
    //        std::unordered_set<uint> e_cuts;
    //        for(uint fid : m.adj_e2f(eid)) if(m.face_data(fid).cut_id>=0) e_cuts.insert(m.face_data(fid).cut_id);
    //        if(( m.edge_is_on_srf(eid) && e_cuts.size()==1) ||
    //           (!m.edge_is_on_srf(eid) && e_cuts.size()==2))
    //        {
    //            ++valence;
    //            v_cuts.insert(e_cuts.begin(), e_cuts.end());
    //        }
    //    }
    //    bool regular     = (valence==0);
    //    bool mm_srf_edge = ( m.vert_is_on_srf(vid) && valence==2 && v_cuts.size()==1) ||
    //                       ( m.vert_is_on_srf(vid) && valence==2 && v_cuts.size()==2); // e.g. two cuts along a concave crease
    //    bool mm_vol_edge = (!m.vert_is_on_srf(vid) && valence==2 && v_cuts.size()==2);
    //    bool mm_srf_vert = ( m.vert_is_on_srf(vid) && valence==5 && v_cuts.size()==2) ||
    //                       ( m.vert_is_on_srf(vid) && valence==4 && v_cuts.size()==3); // e.g. two cuts along a concave crease, plus another cut orthogonal to the crease
    //    bool mm_vol_vert = (!m.vert_is_on_srf(vid) && valence==6 && v_cuts.size()==3);
    //    if(!regular && !mm_srf_edge && !mm_vol_edge && !mm_srf_vert && !mm_vol_vert)
    //    {
    //        std::cout << "Vert " << vid << " is neither regular nor a valid mm_edge or mm_vert. (#incident mm edges: " << valence << ", #cuts: " << v_cuts.size() << ", srf: " << m.vert_is_on_srf(vid) << "). Discard!" << std::endl;
    //        return false;
    //    }
    //}

    return true;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

// verifies that the current cut intersects the other cuts at open curves, and never at a closed loop
//
bool cut_makes_inner_bubbles(const TetMesh                  & m,
                             const std::unordered_set<uint> & edges_on_cut)
{
    std::map<ipair,std::unordered_set<uint>> l_inters; // pair of loop ids => edges in their intersection
    for(uint eid : edges_on_cut)
    {
        // Jan 18 2020: avoid deeming two loops in the same concavity as a bubble
        if(m.edge_is_on_srf(eid)) continue;

        uint count = 0;
        std::unordered_set<uint> cuts;
        for(uint id : m.adj_e2f(eid))
        {
            if(m.face_data(id).cut_id>=0)
            {
                ++count;
                cuts.insert(m.face_data(id).cut_id);
            }
        }
        assert(cuts.size()<=2);
        if(cuts.size()==2) // REMEMBER: bubbles on the srf would not enter this if clause!
        {
            auto  it    = cuts.begin();
            uint  l0    = *it; ++it;
            uint  l1    = *it;
            ipair lpair = unique_pair(l0,l1);
            l_inters[lpair].insert(eid);
        }
    }

    //std::cout << l_inters.size() << " loop intersections found" << std::endl;

    for(auto obj : l_inters)
    {
        std::unordered_set<uint> visited;
        auto & e_inters = obj.second;
        for(uint eid : e_inters)
        {
            if(CONTAINS(visited,eid)) continue;
            visited.insert(eid);
            std::unordered_set<uint> e;
            std::unordered_set<uint> v;
            uint curr = eid;
            while(true)
            {
                e.insert(curr);
                v.insert(m.edge_vert_id(curr,0));
                v.insert(m.edge_vert_id(curr,1));
                bool found = false;
                for(uint nbr : m.adj_e2e(curr))
                {
                    if(CONTAINS(e_inters,nbr) && DOES_NOT_CONTAIN(visited,nbr))
                    {
                        visited.insert(nbr);
                        curr = nbr;
                        found = true;
                        break;
                    }
                }
                if(!found) break;
            }
            if(v.size()==e.size())
            {
                std::cout << "found bubble between loops " << obj.first << std::endl;
                return true;
            }
        }
    }
    return false;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void cut(GlobalState & state, const uint lid)
{
    assert(memory_usage_in_giga_bytes()<10);

    // skip loops already used or previsouly reverted
    if(state.loops.at(lid).used || state.loops.at(lid).reverted || state.loops.at(lid).Nico_bug) return;
    state.loops.at(lid).used = true;

    state.profiler.push("Cut loop #" + std::to_string(lid));

    state.profiler.push("preprocessing");
    std::vector<uint>  vids;        // IDS of loop points
    std::vector<vec3d> points;      // coordinates of loop points
    std::vector<vec3d> srf_normals; // surface normals at each loop point
    std::vector<vec3d> cut_normals; // cut normal at each loop point
    state.loops.trace_loop(lid, vids, points, srf_normals, cut_normals);
    state.profiler.pop("\t" + std::to_string(memory_usage_in_giga_bytes()) + "GB RAM");
    assert(memory_usage_in_giga_bytes()<10);

    state.profiler.push("HRBF");
    Hermite_RBF<CubicRBF> HRBF(points, cut_normals);
    // look for other loops interpolated by the current cut
    if(state.loops.at(lid).type==CONCAVE)
    {
        find_mates(state, HRBF, vids, points, srf_normals, cut_normals);
        HRBF = Hermite_RBF<CubicRBF>(points, cut_normals);
    }
    state.profiler.pop("\t" + std::to_string(memory_usage_in_giga_bytes()) + "GB RAM");
    assert(memory_usage_in_giga_bytes()<10);

    state.profiler.push("Tessellation (" + std::to_string(state.m_vol.num_verts()) + "V)");
    ScalarField f = HRBF.eval(state.m_vol.vector_verts());
    f.copy_to_mesh(state.m_vol);
    for(uint vid : vids) state.m_vol.vert_data(vid).uvw[0] = 0.0; // snap loop verts to level set
    tessellate_connected_component_of_isosurface(state.m_vol, state.loops.at(lid), vids, srf_normals, 0.0);
    state.profiler.pop("\t" + std::to_string(memory_usage_in_giga_bytes()) + "GB RAM");
    assert(memory_usage_in_giga_bytes()<10);

    state.profiler.push("Analyze the cut. Revert if necessary");
    std::unordered_set<uint> faces_on_cut;
    std::unordered_set<uint> edges_on_cut;
    std::unordered_set<uint> verts_on_cut;
    for(uint fid=0; fid<state.m_vol.num_faces(); ++fid)
    {
        uint count = 0;
        for(uint vid : state.m_vol.adj_f2v(fid))
        {
            if(state.m_vol.vert_data(vid).flags[MARKED]) ++count;
        }
        if(count==state.m_vol.adj_f2v(fid).size())
        {
            faces_on_cut.insert(fid);
            for(uint eid : state.m_vol.adj_f2e(fid)) edges_on_cut.insert(eid);
            for(uint vid : state.m_vol.adj_f2v(fid)) verts_on_cut.insert(vid);
        }
    }
    if(!validate_cut_topology(state.m_vol, lid, faces_on_cut, edges_on_cut, verts_on_cut) ||
        cut_makes_inner_bubbles(state.m_vol, edges_on_cut))
    {
        // unmark current cut
        std::cout << "Revert cut" << std::endl;
        revert_cut(state, lid);
        state.profiler.pop();
        state.profiler.pop();
        return;
    }
    state.profiler.pop("\t" + std::to_string(memory_usage_in_giga_bytes()) + "GB RAM");
    assert(memory_usage_in_giga_bytes()<10);

    state.profiler.push("Labeling update");
    uint n = update_poly_labeling(state.m_vol);
    state.profiler.pop("(" + std::to_string(n) + " clusters)");
    if(!batch_mode)
    {
        state.profiler.push("Color wrt label + updateGL");
        state.m_vol.poly_color_wrt_label(false, 0.25, 1.0);
        state.profiler.pop();
    }
    state.profiler.pop("\t" + std::to_string(memory_usage_in_giga_bytes()) + "GB RAM");
    assert(memory_usage_in_giga_bytes()<10);
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void revert_cut(GlobalState & state, const uint lid)
{
    state.loops.at(lid).reverted = true;
    for(uint fid=0; fid<state.m_vol.num_faces(); ++fid)
    {
        if(state.m_vol.face_data(fid).cut_id==(int)lid) state.m_vol.face_data(fid).cut_id = -1;
    }
}
