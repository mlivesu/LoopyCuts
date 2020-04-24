#include "mesh_extractor.h"

MeshExtractor::MeshExtractor(TetMesh & m)
{
    init(m);
    make_verts(m);
    make_edges(m);
    make_faces();
    make_polys(m);

    // just for visualization (to be removed)
    for(uint fid=0; fid<m.num_faces(); ++fid) if(m.face_is_on_srf(fid)) m.face_data(fid).flags[MARKED] = false;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void MeshExtractor::update_v_map()
{
    m2mm_vmap.clear();
    for(uint mm_vid=0; mm_vid<mm.num_verts(); ++mm_vid)
    {
        uint m_vid = mm.vert_data(mm_vid).vid_on_m;
        m2mm_vmap[m_vid] = mm_vid;
    }
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void MeshExtractor::update_f_map()
{
    m2mm_fmap.clear();
    for(uint fid=0; fid<mm.num_faces(); ++fid)
    {
        int l = mm.face_data(fid).label;
        m2mm_fmap[l] = fid;
    }
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void MeshExtractor::init(TetMesh &m)
{
    // mark verts, edges and faces participating in some element of the meta mesh
    // polys of the meta mesh will be clusters of tets having same label, disjoint
    // from the other clusters through sets of such lower dimension simplices
    m.vert_unmark_all();
    m.edge_unmark_all();
    m.face_unmark_all();

    // detect boundary faces (surface faces + inner faces having two incident polys with different labels)
    for(uint fid=0; fid<m.num_faces(); ++fid)
    {
        if(m.face_is_on_srf(fid))
        {
            m.face_data(fid).flags[MARKED] = true;
        }
        else if(m.face_data(fid).cut_id>=0) // the face is on some cut (it doesn't necessarily mean it separates two clusters)
        {
            assert(m.adj_f2p(fid).size()==2);
            uint pid0 = m.adj_f2p(fid).front();
            uint pid1 = m.adj_f2p(fid).back();
            assert(pid0!=pid1);
            if(m.poly_data(pid0).label!=m.poly_data(pid1).label) m.face_data(fid).flags[MARKED] = true;
        }
    }

    // mark edges (either sharp creases or edges having 3+ incident boundary faces)
    for(uint eid=0; eid<m.num_edges(); ++eid)
    {
        if(m.edge_data(eid).must_be_mm_edge())
        {
            m.edge_data(eid).flags[MARKED]=true;
        }
        else
        {
            std::vector<uint> faces;
            for(uint fid : m.adj_e2f(eid))
            {
                if(m.face_data(fid).flags[MARKED]) faces.push_back(fid);
            }
            uint val = faces.size();
            if(val >2) m.edge_data(eid).flags[MARKED] = true;
            // there shouldn't exist any edge having only one incident cut (for surface edges external triangles count as cuts)
            //if(val==1) assert(false && "edges with only one incident cut are not supposed to happen. Something went wrong with the last HRBF cut");
        }
    }

    // mark verts having 1 or 3+ incident edges
    for(uint vid=0; vid<m.num_verts(); ++vid)
    {
        std::vector<uint> edges;
        for(uint eid : m.adj_v2e(vid))
        {
            if(m.edge_data(eid).flags[MARKED]) edges.push_back(eid);
        }
        uint val = edges.size();
        if(val>0 && val!=2) m.vert_data(vid).flags[MARKED] = true;
        // valence one verts may arise in presence of "dead end creases" (e.g. Fandisk lateral "dead end" crease)
        if(val==1) assert(m.edge_data(edges.front()).flags[MARKED] && "valence one verts may happen only at creases");
    }

    // label faces
    m.face_apply_label(-1);
    uint fresh_label = 0;
    for(uint fid=0; fid<m.num_faces(); ++fid)
    {
        if(!m.face_data(fid).flags[MARKED])  continue;
        if(m.face_data(fid).label>=0) continue;
        m.face_data(fid).label = fresh_label;
        bool srf_label = m.face_is_on_srf(fid);
        std::queue<uint> q;
        q.push(fid);
        while(!q.empty())
        {
            uint fid = q.front();
            q.pop();
            for(uint nbr : m.adj_f2f(fid))
            {
                if(!m.face_data(nbr).flags[MARKED]) continue;
                uint eid = m.face_shared_edge(fid,nbr);
                if(m.edge_data(eid).flags[MARKED]) continue;
                if(m.face_data(nbr).label==-1)
                {
                    assert(srf_label==m.face_is_on_srf(nbr));
                    m.face_data(nbr).label = fresh_label;
                    q.push(nbr);
                }
            }
        }
        if(srf_label) srf_face_labels.insert(fresh_label);
        ++fresh_label;
    }
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void MeshExtractor::make_verts(TetMesh & m)
{
    for(uint m_vid=0; m_vid<m.num_verts(); ++m_vid)
    {
        if(m.vert_data(m_vid).flags[MARKED])
        {
            uint mm_vid = mm.vert_add(m.vert(m_vid));
            mm.vert_data(mm_vid).vid_on_m = m_vid;
        }
    }
    update_v_map();
    //std::cout << mm.num_verts() << " MM verts found" << std::endl;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

// WARNING: this will not detect closed edge loops with no verts!
void MeshExtractor::make_edges(TetMesh & m)
{
    std::unordered_set<uint> m_seeds_used;
    uint nv = mm.num_verts();
    for(uint mm_vid=0; mm_vid<nv; ++mm_vid) // nv avoids using the newly inserted vertices as seeds...
    {
        uint m_v_seed = mm.vert_data(mm_vid).vid_on_m;
        for(uint m_e_seed : m.adj_v2e(m_v_seed))
        {
            if(m.edge_data(m_e_seed).flags[MARKED] && DOES_NOT_CONTAIN(m_seeds_used, m_e_seed))
            {
                // extract the ordered chain of vertices participating in the edge
                std::vector<uint> m_v_chain = make_edge(m, m_v_seed, m_e_seed);

                // mark the chain as used at both sides to avoid using it again
                int e_front = m.edge_id(*(m_v_chain.begin()),  *(m_v_chain.begin()+1));  assert(e_front>=0);
                int e_back  = m.edge_id(*(m_v_chain.rbegin()), *(m_v_chain.rbegin()+1)); assert(e_back>=0);
                m_seeds_used.insert(e_front);
                m_seeds_used.insert(e_back);

                uint m_beg  = m_v_chain.front();
                uint m_end  = m_v_chain.back();

                uint mm_beg = m2mm_vmap.at(m_beg);
                uint mm_end = m2mm_vmap.at(m_end);

                auto f_labels = face_labels_incident_to_edge(m, m_e_seed); // faces incident to the edge

                if(m_beg==m_end) continue; // edge is a closed loop: discard
                int mm_eid = mm.edge_id(mm_beg,mm_end);
                if(mm_eid==-1) // new edge
                {
                    mm_eid = mm.edge_add(mm_beg,mm_end);
                    mm.edge_data(mm_eid).v_chain_on_m = m_v_chain; // there may be other ones (for duplicated edges)
                    mm.edge_data(mm_eid).face_labels = f_labels;
                    mm.edge_data(mm_eid).flags[MARKED] = (m.edge_data(m_e_seed).sharp_crease);
                }
                else // edge already existing, just add the new incident face labels
                {
                    for(uint l : f_labels) mm.edge_data(mm_eid).face_labels.insert(l);

                    // if multiple chains connect two mm verts, save all such chains for further analysis
                    auto id = unique_pair(m_beg,m_end);
                    auto it = buggy_chains.find(id);
                    if(it==buggy_chains.end()) buggy_chains[id].push_back(mm.edge_data(mm_eid).v_chain_on_m);
                    buggy_chains.at(id).push_back(m_v_chain);
                }
            }
        }
    }

    // WARNING: this conflicts with the buggy_chains map (vids change...)
    // iteratively remove verts with valence 1 or 2 (and attached edges as well)
    while(remove_degenerate_edges()){}
    // for consistency, in case some verts have been deleted
    update_v_map();

    //std::cout << mm.num_edges() << " MM edges found" << std::endl;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

std::vector<uint> MeshExtractor::make_edge(const TetMesh & m, const uint start_v, const uint start_e) const
{
    uint curr_v = m.vert_opposite_to(start_e, start_v);
    uint curr_e = start_e;
    std::vector<uint> v_chain;
    v_chain.push_back(start_v);
    v_chain.push_back(curr_v);
    while(!m.vert_data(curr_v).flags[MARKED])
    {
        for(uint eid : m.adj_v2e(curr_v))
        {
            if(m.edge_data(eid).flags[MARKED] && eid!=curr_e)
            {
                curr_v = m.vert_opposite_to(eid, curr_v);
                curr_e = eid;
                v_chain.push_back(curr_v);
                break;
            }
        }
    }
    return v_chain;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

bool MeshExtractor::remove_degenerate_edges() // valence 1 verts (and attached edges)
{
    for(uint vid=0; vid<mm.num_verts(); ++vid)
    {
        if(mm.vert_valence(vid)==1)
        {
            //std::cout << mm.num_verts() << " removing valence one vert " << vid << std::endl;
            assert(mm.adj_v2e(vid).size()==1);
            uint eid = mm.adj_v2e(vid).front();
            uint nbr = mm.vert_opposite_to(eid, vid);

            if(nbr>vid) // otherwise deleting vid from mm would change nbr's ID
            {
                mm.vert_switch_id(vid,nbr);
                std::swap(vid,nbr);
            }

            // remove edge from mm and update connectivity
            REMOVE_FROM_VEC(mm.adj_v2v(vid), nbr);
            REMOVE_FROM_VEC(mm.adj_v2v(nbr), vid);
            REMOVE_FROM_VEC(mm.adj_v2e(vid), eid);
            REMOVE_FROM_VEC(mm.adj_v2e(nbr), eid);
            mm.edge_remove_unreferenced(eid);
            mm.vert_remove_unreferenced(vid);
            if(mm.vert_valence(nbr)==0) mm.vert_remove_unreferenced(nbr); else
            if(mm.vert_valence(nbr)==2)
            {
                std::cout << "WARNING: removing a dead end has generated a valence 2 vertex. FixMe" << std::endl;
            }
            return true;
        }
        else if(mm.vert_valence(vid)==2)
        {
            //std::cout << mm.num_verts() << " removing valence two vert " << vid << std::endl;
            assert(mm.adj_v2e(vid).size()==2);

            uint e0 = mm.adj_v2e(vid).front();
            uint e1 = mm.adj_v2e(vid).back();
            if(e0<e1) std::swap(e0,e1);
            uint v0 = mm.vert_opposite_to(e0, vid);
            uint v1 = mm.vert_opposite_to(e1, vid);
            assert(v0!=v1);

            REMOVE_FROM_VEC(mm.adj_v2v(v0),vid);
            REMOVE_FROM_VEC(mm.adj_v2v(v1),vid);
            REMOVE_FROM_VEC(mm.adj_v2e(v0),e0);
            REMOVE_FROM_VEC(mm.adj_v2e(v1),e1);
            if(!mm.verts_are_adjacent(v0,v1))
            {
                uint id = mm.edge_add(v0,v1);
                mm.edge_data(id).face_labels = mm.edge_data(e0).face_labels;
                mm.edge_data(id).flags[MARKED] = mm.edge_data(e0).flags[MARKED];
                for(auto l : mm.edge_data(e1).face_labels) mm.edge_data(id).face_labels.insert(l);
            }
            mm.edge_remove_unreferenced(e0);
            mm.edge_remove_unreferenced(e1);
            mm.vert_remove_unreferenced(vid);
            return true;
        }
    }
    return false;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void MeshExtractor::make_faces()
{
    std::unordered_set<int> visited_labels;
    visited_labels.insert(-1); // default label (to ignore faces not on a cut)

    for(uint mm_eid=0; mm_eid<mm.num_edges(); ++mm_eid)
    {
        for(int l : mm.edge_data(mm_eid).face_labels)
        {
            if(DOES_NOT_CONTAIN(visited_labels, l))
            {
                visited_labels.insert(l);

                std::vector<uint> f = make_face(mm_eid, l);
                if(!f.empty())
                {
                    int id = mm.face_id(f);
                    if(id==-1)
                    {
                        id = mm.face_add(f);
                        mm.face_data(id).label = l;
                    }
                    else if(CONTAINS(srf_face_labels,l))
                    {
                        // in case of duplicates, always keep the reference to the outer face
                        mm.face_data(id).label = l;
                    }
                    m2mm_fmap[l] = (uint)id;
                }
            }
        }
    }

    //std::cout << mm.num_faces() << " MM faces found" << std::endl;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

std::vector<uint> MeshExtractor::make_face(const uint start_e, const int f_label) const
{
    uint curr_e  = start_e;
    uint start_v = mm.edge_vert_id(start_e,0);
    uint curr_v  = mm.edge_vert_id(start_e,1);
    std::vector<uint> f = { start_v };

    do
    {
        f.push_back(curr_v);
        std::vector<uint> next_e;
        for(uint eid : mm.adj_v2e(curr_v))
        {
            if(eid==curr_e) continue;
            if(CONTAINS(mm.edge_data(eid).face_labels,f_label))
            {
                next_e.push_back(eid);
            }
        }
        // face is not sortable: return empty vector
        if(next_e.size()!=1)
        {
            //std::cout << "Unsortable face! (next size: " << next_e.size() << ") - discarded" << std::endl;
            //PRINT(f, "F:");
            //for(uint eid=0; eid<mm.num_edges(); ++eid)
            //{
            //    PRINT(mm.edge_data(eid).face_labels, "flabels for " + std::to_string(eid));
            //}
            return std::vector<uint>();
        }
        curr_e = next_e.front();
        curr_v = mm.vert_opposite_to(curr_e,curr_v);
    }
    while(curr_v!=start_v && f.size()<mm.num_verts()); // to avoid inf loops
    assert(curr_v==start_v);
    return f;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

std::unordered_set<uint> MeshExtractor::face_labels_incident_to_edge(const TetMesh & m, const uint eid)
{
    std::unordered_set<uint> labels;
    for(uint fid : m.adj_e2f(eid))
    {
        if(m.face_data(fid).flags[MARKED])
        {
            int l = m.face_data(fid).label;
            assert(l>=0);
            labels.insert(l);
        }
    }
    return labels;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void MeshExtractor::make_polys(TetMesh & m)
{
    uint discarded = 0;
    m.poly_unmark_all();
    for(uint pid=0; pid<m.num_polys(); ++pid)
    {
        if(!m.poly_data(pid).flags[MARKED])
        {
            std::unordered_set<uint> p = make_poly(m, pid);
            std::vector<uint>        f(p.begin(), p.end());

            if(topological_checks(p) && mm.poly_id(f)==-1)
            {
                std::vector<bool> w(p.size(), true);
                mm.poly_add(f,w);
            }
            else ++discarded;
        }
    }

    //std::cout << mm.num_polys() << " MM polys found (" << discarded << " clusters were discarded)" << std::endl;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

std::unordered_set<uint> MeshExtractor::make_poly(TetMesh & m, const uint start_p) const
{
    assert(!m.poly_data(start_p).flags[MARKED]);
    m.poly_data(start_p).flags[MARKED] = true;
    std::queue<uint> q;
    q.push(start_p);
    std::vector<bool> visited(m.num_faces(), false);
    std::unordered_set<uint> flist;
    while(!q.empty())
    {
        uint pid = q.front();
        q.pop();
        for(uint fid : m.adj_p2f(pid))
        {
            if(visited.at(fid)) continue;
            visited.at(fid) = true;

            int l = m.face_data(fid).label;
            if(m.face_data(fid).flags[MARKED])
            {
                auto mm_fid_it = m2mm_fmap.find(l);
                // need to check whether the label has a meta mesh face id or not
                // (faces with no incident vertices and dangling faces are discarded by default)
                if(mm_fid_it!=m2mm_fmap.end())
                {
                    flist.insert(mm_fid_it->second);
                }
            }
            else
            {
                int nbr = m.poly_adj_through_face(pid, fid);
                //assert(nbr>=0);
                if(nbr>=0 && !m.poly_data(nbr).flags[MARKED])
                {
                    m.poly_data(nbr).flags[MARKED] = true;
                    q.push(nbr);
                }
            }
        }
    }
    return flist;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

bool MeshExtractor::topological_checks(const std::unordered_set<uint> & p)
{
    if(p.empty()) return false;

    // collect all vertices and edges
    std::unordered_set<uint> v;
    std::unordered_set<uint> e;
    for(uint fid : p)
    for(uint eid : mm.adj_f2e(fid))
    {
        e.insert(eid);
        v.insert(mm.edge_vert_id(eid,0));
        v.insert(mm.edge_vert_id(eid,1));
    }

    // WATERTIGHTNESS CHECK: each edge should have exactly 2 incident faces
    for(uint eid : e)
    {
        std::vector<uint> f;
        for(uint fid : mm.adj_e2f(eid))
        {
            if(CONTAINS(p, fid)) f.push_back(fid);
        }
        if(f.size()!=2)
        {
            //++polys_non_watertight;
            //std::cout << "WARNING: non watertight poly!" << std::endl;
            //PRINT(f, "faces incident to edge " + std::to_string(eid));
            return false;
        }
    }

    // MANIFOLDNESS CHECKS: each face should have either 1 or 2 incident polys
    for(uint fid : p)
    {
        if(mm.adj_f2p(fid).size()>=2)
        {
            ++polys_non_manifold;
            //std::cout << "WARNING: face has already 2 incident polys!" << std::endl;
            //PRINT(mm.adj_f2p(fid), "polys incident to face " + std::to_string(fid));
            return false;
        }
    }

    // GENUS CHECK: each poly should be a topological ball (genus>0 is incompatible with midpoint)
    if(v.size() - e.size() + p.size() != 2)
    {
        ++polys_with_high_genus;
        //std::cout << "WARNING: poly has not genus zero!" << std::endl;
        return false;
    }

    return true;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

bool MeshExtractor::converged(const TetMesh & m) const
{
    if(topological_convergence(m) && geometric_convergence(m))
    {
        std::cout << "CONVERGED!" << std::endl;
        return true;
    }
    return false;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

bool MeshExtractor::topological_convergence(const TetMesh & m) const
{
    int g0 = mm.genus();
    int g1 = m.genus();

    //std::cout << "::::::::::::: TOPOLOGICAL CHECKS :::::::::::::" << std::endl;
    //std::cout << g0                      << " genus (tetmesh genus is: " << g1 << ")          " << std::endl;
    //std::cout << polys_non_manifold      << " polys non manifold (face with 3+ incident polys)" << std::endl;
    //std::cout << polys_with_high_genus   << " polys with high genus                           " << std::endl;
    //std::cout << "::::::::::::::::::::::::::::::::::::::::::::::" << std::endl;

    if(g0==g1 && polys_non_manifold==0 && polys_with_high_genus==0) return true;
    return false;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

bool MeshExtractor::geometric_convergence(const TetMesh & m) const
{
    double worst_ratio = 1.0;
    for(uint eid=0; eid<mm.num_edges(); ++eid)
    {
        if(mm.adj_e2p(eid).empty()) continue;
        double edge_length  = mm.edge_length(eid);
        double chain_length = 0;
        auto & chain = mm.edge_data(eid).v_chain_on_m;
        if(chain.empty()) continue;
        for(uint i=0; i<chain.size()-1; ++i)
        {
            chain_length += m.vert(chain.at(i)).dist(m.vert(chain.at(i+1)));
        }
        worst_ratio = std::min(worst_ratio, edge_length/chain_length);
    }

    double target_area = m.mesh_srf_area();
    double area        = mm.mesh_srf_area();
    double coverage    = std::fabs(target_area-area)/target_area;

    // to avoid selecting orrible midpoints in midpoint refinement,
    // make sure non quad faces have a good aspect ratio
    for(uint fid=0; fid<mm.num_faces(); ++fid)
    {
        if(mm.face_is_on_srf(fid) && mm.verts_per_face(fid)!=4)
        {
            double min_l = inf_double;
            double max_l = 0;
            for(uint eid : mm.adj_f2e(fid))
            {
                bool v0_locked = false;
                bool v1_locked = false;
                uint v0 = mm.edge_vert_id(eid,0);
                uint v1 = mm.edge_vert_id(eid,1);
                for(uint e : mm.adj_v2e(v0)) if(mm.edge_data(e).flags[MARKED]) { v0_locked = true; break; }
                for(uint e : mm.adj_v2e(v1)) if(mm.edge_data(e).flags[MARKED]) { v1_locked = true; break; }
                // restrict ratio chck only to edges where both endpoints participate in a crease and are locked.
                // otherwise trust the smoother to let vertices arrange in a way that gives mesh regularity
                if(v0_locked && v1_locked)
                {
                    auto l = mm.edge_length(eid);
                    min_l = std::min(min_l, l);
                    max_l = std::max(max_l, l);
                }
            }
            if(max_l>0 && min_l<inf_double && max_l/min_l>10) return false;
        }
    }

    //std::cout << "::::::::::::: GEOMETRIC CHECKS :::::::::::::::" << std::endl;
    //std::cout << target_area << " tetmesh area                  " << std::endl;
    //std::cout << area        << " meta mesh area                " << std::endl;
    //std::cout << coverage    << " coverage                      " << std::endl;
    //std::cout << worst_ratio << " worst edge ratio              " << std::endl;
    //std::cout << "::::::::::::::::::::::::::::::::::::::::::::::" << std::endl;

    if(coverage<0.1) return true;
    return false;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void MeshExtractor::post_convergence_cut_analysis(const TetMesh            & m,
                                                  const Loops              & /*loops*/,
                                                  std::unordered_set<uint> & cuts_to_do,     // cuts that would fix some existing val2 verts
                                                  std::unordered_set<uint> & cuts_to_undo)  // cuts that created val2 verts that cannot be removed
{
    std::cout << "\nPost Convergence Cut Analysis" << std::endl;

    cuts_to_do.clear();
    cuts_to_undo.clear();

    // go through each MM edge connected by more than one chain
    for(auto e : buggy_chains)
    {
        /* for each such defect, consider all the chains of edges in M that
         * connect two endpoints of a MM edge. A subset of these chains will
         * be formed by surface edges of M.
         *
         *  - take all these chains and find the set of loops that intersect them all.
         *    Cutting along one of these loops should fix the bug.
         *
         *  - if no crossing loops are found, the only way to fix is to revert all
         *    cuts that created such chains but one
         *
         *  - finally, if there are no surface chain, it means this is sort of an internal
         *    pocket made by 3 or more cuts that intersect in a bad way. the solution is:
         *    revert at least one of such cuts
        */

        if(DOES_NOT_CONTAIN(m2mm_vmap, e.first.first)) continue;
        if(DOES_NOT_CONTAIN(m2mm_vmap, e.first.second)) continue;
        uint vid0    = m2mm_vmap.at(e.first.first);
        uint vid1    = m2mm_vmap.at(e.first.second);
        int  mm_eid  = mm.edge_id(vid0,vid1);
        bool srf_bug = mm.vert_is_on_srf(vid0) && mm.vert_is_on_srf(vid1);
        std::cout << " MM edge (SRF:" << srf_bug << ") "<< mm_eid << e.first << " is connected through " << e.second.size() << " chains" << std::endl;
        assert(mm_eid>=0);

        // for each chain, report whether it is on the surface, and what loops it intersects
        bool has_surface_chains = false;
        std::vector<bool> chain_is_on_srf;
        std::vector<std::vector<uint>> chain_crosses_loops;
        auto & chains = e.second;
        for(auto chain : chains)
        {
            assert(chains.size()>1);
            int m_eid = m.edge_id(chain.at(0), chain.at(1));
            assert(m_eid>=0);

            if(!m.edge_is_on_srf(m_eid))
            {
                std::cout << "\tinner chain" << std::endl;
                chain_is_on_srf.push_back(false);
                chain_crosses_loops.push_back(std::vector<uint>());
            }
            else
            {
                std::vector<uint> loops_crossed;                
                int lid = m.edge_data(m_eid).loop_ids[0]; // this may also be -1 (in case the chain comes from a loop that was not in the original loop set)
                for(uint i=1; i<chain.size()-1; ++i)
                {
                    uint vid = chain.at(i);
                    for(uint eid : m.vert_adj_srf_edges(vid))
                    {
                        int id = m.edge_data(eid).loop_ids[0];
                        if(id>=0 && id!=lid)
                        {
                            loops_crossed.push_back(id);
                            break;
                        }
                    }
                }
                has_surface_chains = true;
                chain_is_on_srf.push_back(true);
                chain_crosses_loops.push_back(loops_crossed);
                PRINT(loops_crossed, "\tsurface chain, intersecting loops");
            }
        }

        // Surface issue. Find the loops that cross all surface chains, and cut along them to fix it
        if(has_surface_chains)
        {
            // intersect all loops to find one that traverses all chains
            bool found_first = false;
            std::vector<uint> res;
            for(uint i=0; i<chain_is_on_srf.size(); ++i)
            {
                if(chain_is_on_srf.at(i))
                {
                    if(!found_first)
                    {
                        found_first = true;
                        res = chain_crosses_loops.at(i);
                    }
                    else
                    {
                        auto tmp = res;
                        SET_INTERSECTION(chain_crosses_loops.at(i), tmp, res, true);
                    }
                }
            }
            if(!res.empty())
            {
                PRINT(res, "Surface Issue. Can fix cutting along loops");
                //cuts_to_do.insert(res.begin(), res.end());
                cuts_to_do.insert(*res.begin());
            }
            else // if there is no such a loop, revert a cut to fix things
            {
                std::cout << "Surface Issue. Couldn't find loops to cut through to fix things. Revert one cut" << std::endl;
                int eid = m.edge_id(chains.front().at(0), chains.front().at(1));
                for(uint fid : m.adj_e2f(eid))
                {
                    if(m.face_data(fid).cut_id>=0)
                    {
                        cuts_to_undo.insert(m.face_data(fid).cut_id);
                        break;
                    }
                }
            }
        }
        else
        {
            std::cout << "Inner issue. Revert one cut" << std::endl;
            int eid = m.edge_id(chains.front().at(0), chains.front().at(1));
            for(uint fid : m.adj_e2f(eid))
            {
                if(m.face_data(fid).cut_id>=0)
                {
                    cuts_to_undo.insert(m.face_data(fid).cut_id);
                    break;
                }
            }
        }
    }
}

