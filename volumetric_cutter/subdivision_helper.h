#ifndef SUBDIVISION_HELPER_H
#define SUBDIVISION_HELPER_H

#include "definitions.h"
#include <unordered_map>
#include <cinolib/octree.h>


//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

typedef struct
{
    Trimesh<> m;                              // mesh of the local patch
    std::vector<uint>  corners;               // patch corners (local vids)
    std::unordered_map<uint,uint> mm2corners; // connects MM verts with local corners
    bool flawed = false;                      // if anything unexpected occurs... discard it
    Octree octree = Octree(4,100);
}
Patch;

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

class SubdivisionHelper
{
    public:

        SubdivisionHelper(const TetMesh  & m,
                          const MetaMesh & mm);

        //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        MetaMesh subdivide();

        //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        void extract_surface_patches();
        void make_uv_maps();
        bool map_to_polygon(Trimesh<> & m, const std::vector<uint> & corners);

        //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    private:

        MetaMesh mm_refined;
        const TetMesh  & m;
        const MetaMesh & mm;

        std::unordered_map<int,Patch> patches; // one for each MM face
};

#endif // SUBDIVISION_HELPER_H
