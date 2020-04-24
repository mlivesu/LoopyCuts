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
