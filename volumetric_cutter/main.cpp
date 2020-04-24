#include <QApplication>
#include <cinolib/gui/qt/qt_gui_tools.h>
#include <cinolib/tetgen_wrap.h>
#include "definitions.h"
#include "state.h"
#include "batch.h"
#include "gui.h"

using namespace cinolib;

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

bool        batch_mode;
std::string model_name;
std::string loops_name;
std::string output_name;

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void parse_command_line(int argc, char **argv)
{
    assert(argc>=3);
    model_name = std::string(argv[1]);
    loops_name = std::string(argv[2]);
    batch_mode = (argc==5 && strcmp(argv[3],"-batch-mode")==0);
    if(batch_mode) output_name = argv[4];

    std::cout << std::endl;
    std::cout << "input mesh  : " << model_name << std::endl;
    std::cout << "input loops : " << loops_name << std::endl;
    std::cout << "batch mode  : " << batch_mode << std::endl;
    if(batch_mode)
    {
        std::cout << "output mesh : " << output_name << std::endl;
    }
    std::cout << std::endl;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

TetMesh tetrahedralize(const SrfMesh & m, const double target_num_elements)
{
    std::vector<double> coords;
    std::vector<uint>   edges, tets;
    double max_tet_vol = m.mesh_volume()/target_num_elements;  //std::pow(target_num_elements,3.0)/(6.0*sqrt(2.0));
    // https://stackoverflow.com/questions/16605967/set-precision-of-stdto-string-when-converting-floating-point-values
    if(max_tet_vol<1e-5) max_tet_vol = 1e-5;
    std::string flags  = "QYqa" + std::to_string(max_tet_vol);
    //std::cout << flags << std::endl;
    tetgen_wrap(serialized_xyz_from_vec3d(m.vector_verts()), serialized_vids_from_polys(m.vector_polys()), edges, flags.c_str(), coords, tets);
    return TetMesh(coords, tets);
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void refine_mesh(TetMesh & m)
{
    uint e_count = 0;
    uint f_count = 0;
    uint p_count = 0;

    std::vector<ipair> edges_to_split;
    for(uint eid=0; eid<m.num_edges(); ++eid)
    {
        uint vid0 = m.edge_vert_id(eid,0);
        uint vid1 = m.edge_vert_id(eid,1);

        if(m.edge_data(eid).loop_ids[0]==-1 &&
           m.vert_data(vid0).on_loop        &&
           m.vert_data(vid1).on_loop)
        {
            edges_to_split.push_back(std::make_pair(m.edge_vert_id(eid,0),
                                                    m.edge_vert_id(eid,1)));
        }
    }
    for(auto e : edges_to_split)
    {
        int eid = m.edge_id(e);
        m.edge_split(eid);
        ++e_count;
    }

    std::vector<std::vector<uint>> faces_to_split;
    for(uint fid=0; fid<m.num_faces(); ++fid)
    {
        uint vid0 = m.face_vert_id(fid,0);
        uint vid1 = m.face_vert_id(fid,1);
        uint vid2 = m.face_vert_id(fid,2);

        if(m.vert_data(vid0).on_loop &&
           m.vert_data(vid1).on_loop &&
           m.vert_data(vid2).on_loop)
        {
            faces_to_split.push_back(m.face_verts_id(fid));
        }
    }
    for(auto f : faces_to_split)
    {
        int fid = m.face_id(f);
        m.face_split(fid);
        ++f_count;
    }

    std::vector<std::vector<uint>> tets_to_split;
    for(uint pid=0; pid<m.num_polys(); ++pid)
    {
        uint vid0 = m.poly_vert_id(pid,0);
        uint vid1 = m.poly_vert_id(pid,1);
        uint vid2 = m.poly_vert_id(pid,2);
        uint vid3 = m.poly_vert_id(pid,3);

        if(m.vert_data(vid0).on_loop &&
           m.vert_data(vid1).on_loop &&
           m.vert_data(vid2).on_loop &&
           m.vert_data(vid3).on_loop)
        {
            tets_to_split.push_back(m.poly_faces_id(pid));
        }
    }
    for(auto p : tets_to_split)
    {
        int pid = m.poly_id(p);
        m.poly_split(pid);
        ++p_count;
    }

    std::cout << "Refined mesh to fully detach cuts" << std::endl;
    std::cout << e_count << " edges were split" << std::endl;
    std::cout << f_count << " faces were split" << std::endl;
    std::cout << p_count << " polys were split" << std::endl;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

int main(int argc, char **argv)
{
    if(argc!=3 && argc!=5)
    {
        std::cout << "                                                           " << std::endl;
        std::cout << "Usage:                                                     " << std::endl;
        std::cout << "                                                           " << std::endl;
        std::cout << "  ./volumetric_cut <mesh> <loops> [-batch-mode <folder>]   " << std::endl;
        std::cout << "                                                           " << std::endl;
        std::cout << "  mesh   is a triangle mesh, either in OFF or OBJ format   " << std::endl;
        std::cout << "  loops  is a text file (see the examples for the format)  " << std::endl;
        std::cout << "                                                           " << std::endl;
        std::cout << "  NOTE: if batch-mode is activated, the program runs with  " << std::endl;
        std::cout << "  no GUI and processes ALL the loops in the given order,   " << std::endl;
        std::cout << "  until convergence. A polyhedral mesh will be produced in " << std::endl;
        std::cout << "  output and saved in <folder>.                            " << std::endl;
        std::cout << "                                                           " << std::endl;
        exit(0);
    }
    parse_command_line(argc,argv);

    GlobalState state; // this contains ALL the application data, gui, profiler...
    state.m_srf = SrfMesh(model_name.c_str());
    state.m_vol = tetrahedralize(state.m_srf, 20000);
    state.loops = Loops(loops_name.c_str(), &state.m_srf, &state.m_vol);    
    refine_mesh(state.m_vol);
    state.m_vol.poly_apply_labels(std::vector<int>(state.m_vol.num_polys(),0));

    if(batch_mode)
    {
        run_batch(state);
        return 0;
    }
    else
    {
        state.history.push_state(State(state.m_vol, state.loops));

        QApplication app(argc, argv);
        GUI gui;
        init_gui(gui, state);
        init_events(gui, state);
        // CMD+1 to show trimesh  controls.
        // CMD+2 to show tetmesh  controls.
        // CMD+3 to show polymesh controls.
        SurfaceMeshControlPanel<SrfMesh> tri_panel(&state.m_srf,   &gui.left);
        VolumeMeshControlPanel<TetMesh>  tet_panel(&state.m_vol,   &gui.left);
        VolumeMeshControlPanel<MetaMesh> poly_panel(&state.m_poly, &gui.right);
        QApplication::connect(new QShortcut(QKeySequence(Qt::CTRL+Qt::Key_1), &gui.window), &QShortcut::activated, [&](){tri_panel.show();});
        QApplication::connect(new QShortcut(QKeySequence(Qt::CTRL+Qt::Key_2), &gui.window), &QShortcut::activated, [&](){tet_panel.show();});
        QApplication::connect(new QShortcut(QKeySequence(Qt::CTRL+Qt::Key_3), &gui.window), &QShortcut::activated, [&](){poly_panel.show();});
        return app.exec();
    }
}
