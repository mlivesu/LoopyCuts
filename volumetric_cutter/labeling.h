#ifndef LABELING_H
#define LABELING_H

#include "definitions.h"

uint update_poly_labeling(TetMesh & m)
{
    m.poly_apply_label(-1);
    uint fresh_label = 0;
    for(uint pid=0; pid<m.num_polys(); ++pid)
    {
        if(m.poly_data(pid).label>=0) continue;
        m.poly_data(pid).label = fresh_label;
        std::queue<uint> q;
        q.push(pid);
        while(!q.empty())
        {
            uint pid = q.front();
            q.pop();
            for(uint fid : m.adj_p2f(pid))
            {
                if(m.face_data(fid).cut_id>=0) continue;
                int nbr = m.poly_adj_through_face(pid, fid);
                if(nbr>=0 && m.poly_data(nbr).label==-1)
                {
                    m.poly_data(nbr).label = fresh_label;
                    q.push(nbr);
                }
            }
        }
        ++fresh_label;
    }
    return fresh_label;
}

#endif // LABELING_H
