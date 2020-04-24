#include "loops.h"
#include <fstream>
#include <cinolib/min_max_inf.h>
#include <cinolib/stl_container_utilities.h>
#include <cinolib/gl/draw_arrow.h>

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

Loops::Loops(const char *filename, SrfMesh *m_srf, TetMesh *m_vol)
    : m_srf(m_srf)
    , m_vol(m_vol)
{
    if(m_srf == NULL) return;

    length = m_srf->edge_avg_length()*2.0;
    thick  = length * 0.1;

    FILE *f = fopen(filename, "r");
    assert(f!=nullptr);

    uint num_loops;
    assert(eat_uint(f, num_loops));
    std::cout << num_loops << " loops found" << std::endl;

    m_srf->edge_unmark_all();

    for(uint lid=0; lid<num_loops; ++lid)
    {
        Loop l;

        char word[32];
        assert(eat_word(f, word));
        if(strcmp(word,"REGULAR")==0) l.type = REGULAR; else
        if(strcmp(word,"CONCAVE")==0) l.type = CONCAVE; else
        if(strcmp(word,"CONVEX")==0)  l.type = CONVEX;  else
        assert(false);

        assert(eat_word(f, word));
        if(strcmp(word,"Closed")==0) l.closed = true;  else
        if(strcmp(word,"Open")==0)   l.closed = false; else
        assert(false);

        // a flawed loop is a loop that has two consecutive intersections with
        // some other loop. Cutting along it will produce non hexa elements.
        // we still need to figure out whether we want to use them or not....
        assert(seek_keyword(f, "Cross"));
        assert(eat_word(f, word));
        if(strcmp(word,"OK")==0)   l.flawed = false; else
        if(strcmp(word,"FAIL")==0) l.flawed = true;  else
        assert(false);

        uint num_segs;
        assert(eat_uint(f, num_segs));

        std::cout << "Loop " << lid << "\t" << num_segs << " segs\t" << TYPES_STR[l.type] << "\t" << ((l.closed)?"Closed\t":"Open\t") << ((l.flawed)?"Flawed":"") << std::endl;

        std::map<uint,uint> vert_valence;
        for(uint i=0; i<num_segs; ++i)
        {
            // read pid, edge local ID and bool flag for loop normal (single/double sided)
            uint pid, offset, one_sided_normal;
            assert(eat_uint(f, pid));
            assert(eat_uint(f, offset));
            assert(eat_uint(f, one_sided_normal));

            if(l.Nico_bug) continue; // I am already discarding this loop, I just need to read it all from file....

            uint srf_eid = m_srf->poly_edge_id(pid, offset);
            uint vid0    = m_srf->edge_vert_id(srf_eid,0);
            uint vid1    = m_srf->edge_vert_id(srf_eid,1);
             int vol_eid = m_vol->edge_id(vid0,vid1);
            assert(vol_eid>=0);

            if(l.type==CONVEX && m_vol->edge_data(vol_eid).type==CONVEX) // NICO's BUG: two convex features passing through the same edge (duplicated loops)
            {
                std::cout << "WARNING: I am skipping loop " << lid << " because it seems to be a copy of loop " << m_vol->edge_data(vol_eid).loop_ids[0] << std::endl;
                l.Nico_bug = true;
                continue;
            }

            m_vol->edge_data(vol_eid).type = l.type;
            m_vol->vert_data(vid0).on_loop = true;
            m_vol->vert_data(vid1).on_loop = true;

            // handle loop srf normals
            vec3d n = m_srf->poly_data(pid).normal;
            if(one_sided_normal==0) // average tri normals at both sides of eid
            {
                int opp = m_srf->poly_opposite_to(srf_eid,pid);
                assert(opp>=0);
                n += m_srf->poly_data(opp).normal;
                n.normalize();
                l.polys.insert(opp); // for manual picking
            }
            l.polys.insert(pid); // for manual picking

            if(m_vol->edge_data(vol_eid).loop_ids[0]==-1) // regular (one loop insisting on the edge)
            {
                m_vol->edge_data(vol_eid).loop_ids[0]=lid;
                m_vol->edge_data(vol_eid).srf_normals[0] = n;
                if(l.type==CONVEX)
                {
                    m_vol->edge_data(vol_eid).sharp_crease = true;
                    m_srf->edge_data(srf_eid).flags[MARKED] = true; // used by the smoother
                }
            }
            else // concavity (edge shared by two loops)
            {
                assert(m_vol->edge_data(vol_eid).loop_ids[1]==-1);
                m_vol->edge_data(vol_eid).loop_ids[1]=lid;
                m_vol->edge_data(vol_eid).srf_normals[1] = n;
                m_vol->edge_data(vol_eid).sharp_crease = true;
                m_srf->edge_data(srf_eid).flags[MARKED] = true; // used by the smoother
            }

            // per loop vert valences (used to find startng point for open loops)
            for(uint vid : m_srf->adj_e2v(srf_eid))
            {
                if(CONTAINS(vert_valence,vid))
                {
                    assert(vert_valence.at(vid)==1 && "Loop self intersects!");
                    vert_valence.at(vid)++;
                }
                else vert_valence[vid] = 1;
            }

            // rendering info
            l.drawlist.push_back(m_vol->vert(vid0));
            l.drawlist.push_back(m_vol->vert(vid1));
        }

        if(!l.Nico_bug)
        {
            // find a seed for loop tracing (for open curves always start from a valence 1 vert)
            l.starting_vert = vert_valence.begin()->first;
            if(!l.closed)
            {
                for(auto it : vert_valence)
                {
                    if(it.second==1)
                    {
                        l.starting_vert = it.first;
                        break;
                    }
                }
            }
        }
        // add the new loop
        push_back(l);
    }
    fclose(f);
    assert(size()==num_loops);

    close_open_creases();
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void Loops::draw(const float) const
{
    glDisable(GL_LIGHTING);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_BLEND);
    glDepthRange(0.0, 1.0);
    glDepthFunc(GL_LEQUAL);

    for(uint lid=0; lid<this->size(); ++lid)
    {
        const Loop & l = this->at(lid);

        if(l.active)             glColor3fv(Color::RED().rgba);      else
        if(l.type==CONCAVE)      glColor3fv(Color::BLUE().rgba);     else
        if(l.type==REGULAR)      glColor3fv(Color::GRAY().rgba);     else
        if(l.type==CONVEX)       glColor3fv(Color::GREEN(.8).rgba);  else
        if(l.type==TOP_RELEVANT) glColor3fv(Color::YELLOW(.8).rgba); else
        assert(false);

        glLineWidth(l.active ? 10.f : 2.f);

        if(show_loops || (l.active && show_active_loop))
        {
            for(auto i=l.drawlist.begin(),j=i+1; i<l.drawlist.end(); i+=2,j+=2)
            {
                glBegin(GL_LINES);
                glVertex3dv(i->ptr());
                glVertex3dv(j->ptr());
                glEnd();
            }
        }

        if(l.active && show_active_loop_normals)
        {
            std::vector<uint>  vids;
            std::vector<vec3d> points;
            std::vector<vec3d> srf_normals;
            std::vector<vec3d> cut_normals;
            trace_loop(lid, vids, points, srf_normals, cut_normals);

            // SRF NORMALS
            for(uint i=0; i<points.size(); ++i)
            {
                vec3d A = points.at(i);
                vec3d B = A + srf_normals.at(i) * length;
                arrow(A, B, thick, Color::GREEN().rgba);
            }
            //// CUT TANGENT
            //for(uint i=0; i<points.size(); ++i)
            //{
            //    vec3d A = points.at(i);
            //    vec3d B = points.at((i+1)%points.size());
            //    vec3d C = points.at((i+2)%points.size());
            //    vec3d t = C -A;
            //    t.normalize();
            //    t*= length;
            //    arrow(B, B+t, thick, Color::RED().rgba);
            //}
            // CUT NORMALS
            for(uint i=0; i<points.size(); ++i)
            {
                vec3d A = points.at(i);
                vec3d B = A + cut_normals.at(i) * length;
                arrow(A, B, thick, Color::BLUE().rgba);
            }
        }
    }
    glDepthFunc(GL_LESS);
    glEnable(GL_LIGHTING);
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void Loops::set_active_loop(const uint curr, const bool deactivate_others)
{
    if(deactivate_others) deactivate_all();
    this->at(curr).active = true;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void Loops::activate_all()
{
    for(uint i=0; i<size(); ++i)
    {
        this->at(i).active = true;
    }
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void Loops::deactivate_all()
{
    for(uint i=0; i<size(); ++i)
    {
        this->at(i).active = false;
    }
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

uint Loops::pick_loop(const vec3d & p) const
{
    uint   best = 0;
    double dist = inf_double;
    for(uint lid=0; lid<size(); ++lid)
    {
        double tmp = inf_double;
        for(uint pid : this->at(lid).polys)
        {
            tmp = std::min(tmp, m_srf->poly_centroid(pid).dist(p));
            if(tmp < dist)
            {
                dist = tmp;
                best = lid;// m_srf->poly_data(pid).loop_to_pick;
            }
        }
    }
    return best;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void Loops::trace_loop(const uint           lid,
                       std::vector<uint>  & vids,
                       std::vector<vec3d> & points,
                       std::vector<vec3d> & srf_normals,
                       std::vector<vec3d> & cut_normals) const
{
    uint starting_vert = this->at(lid).starting_vert;

    vids.clear();
    vids.push_back(starting_vert);
    points.clear();
    points.push_back(m_vol->vert(starting_vert));
    srf_normals.clear();
    cut_normals.clear();

    int prev = -1;
    while(true)
    {
        int eid = -1;
        for(uint id : m_vol->adj_v2e(vids.back()))
        {
            if(!m_vol->edge_data(id).belongs_to_loop(lid) || prev==(int)id) continue;
            eid = id;
            break;
        }

        if(eid>=0)
        {
            uint vid = m_vol->vert_opposite_to(eid, vids.back());
            if(vid==starting_vert)
            {
                assert(this->at(lid).closed);
                srf_normals.push_back(srf_normals.back());
                cut_normals.push_back(cut_normals.back());
                break; // closed loop
            }
            vec3d p = m_vol->vert(vid);
            vec3d n = m_vol->edge_data(eid).srf_normal(lid);
            vec3d t = p - points.back();
            vec3d c = t.cross(n);
            c.normalize();
            vids.push_back(vid);
            points.push_back(p);
            srf_normals.push_back(n);
            cut_normals.push_back(c);
            prev = eid;
        }
        else // open loop
        {
            assert(!this->at(lid).closed);
            srf_normals.push_back(srf_normals.back());
            cut_normals.push_back(cut_normals.back());
            break;
        }
    }

    assert(points.size() == cut_normals.size());
    assert(points.size() == srf_normals.size());
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void Loops::close_open_creases()
{
    // Mark as crease any loop that closes a convex crease with a T-junction
    for(uint vid=0; vid<m_vol->num_verts(); ++vid)
    {
        if(!m_vol->vert_is_on_srf(vid)) continue;

        uint count = 0;
        uint lid;
        for(uint eid : m_vol->vert_adj_srf_edges(vid))
        {
            if(m_vol->edge_data(eid).type==CONVEX || m_vol->edge_data(eid).type==CONCAVE)
            {
                ++count;
                //assert(!m_vol->edge_data(eid).has_two_loops()); // false for concave loops
                lid = m_vol->edge_data(eid).loop_ids[0];
            }
        }

        /* Dead end convex crease detected. To avoid making any non removable valence 2 vert in the
         * meta mesh, find the regular loop that "closes the lid", and make sure it will be included
         * in the meta mesh as surface edge, without being cut through
         *
         * WARNING: when such a "lid" loop intersects a concave loop that is cut through twice,
         * the blocks that see pieces of the concave loop as a surface edge (without having surface
         * faces incident to it) will see the intersection point between the lid and the concave loop
         * as a valence 2. A proper splitting scheme must be devised in order to avoid such pathological cases
        */
        if(count==1)
        {
            std::cout << "Vertex " << vid << " is the endpoint of open sharp crease line #" << lid << std::endl;

            std::unordered_set<uint> lid_loops;
            for(uint eid : m_vol->vert_adj_srf_edges(vid))
            {
                int l0 = m_vol->edge_data(eid).loop_ids[0];
                int l1 = m_vol->edge_data(eid).loop_ids[1];
                if(l0!=(int)lid && l0!=-1) lid_loops.insert(l0);
                if(l1!=(int)lid && l1!=-1) lid_loops.insert(l1);
            }
            PRINT(lid_loops, "\tLid loops");
            // if there are two cuts that close the dead end everythng is fine
            // (as long as I can actudally cut through both loops) [example: joint, fandisk]
            //
            // UPDATE 10 Dec 2019: fandisk is a clear example that I need to mark as topologically
            // relevant also concave loops passing from there! the reason is that if I am able to
            // cut twice along a concavity then everything is fine. But if I fail one cut, then I
            // need the second to exist at least on the surface of the meat mesh, to obtain valence 3 everywhere...
            //if(lid_loops.size()==1)
            for(uint lid : lid_loops)
            {
                //int lid = *lid_loops.begin();
                Loop & l = this->at(lid);
                l.used = true; // makes sure no one will cut through it
                l.type = TOP_RELEVANT;

                std::vector<uint>  vids;        // IDS of loop points
                std::vector<vec3d> points;      // coordinates of loop points
                std::vector<vec3d> srf_normals; // surface normals at each loop point
                std::vector<vec3d> cut_normals; // cut normal at each loop point
                trace_loop(lid, vids, points, srf_normals, cut_normals);

                for(auto i=vids.begin(),j=i+1; j<vids.end(); ++i,++j)
                {
                    int eid = m_vol->edge_id(*i,*j);
                    assert(eid>=0);
                    m_vol->edge_data(eid).topologically_needed = true; // makes sure these edges will form edges in the meta mesh
                }
                if(l.closed)
                {
                    int eid = m_vol->edge_id(*vids.begin(),*vids.rbegin());
                    assert(eid>=0);
                    m_vol->edge_data(eid).topologically_needed = true;
                }
            }
        }
    }
}
