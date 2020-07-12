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

#include "gui.h"
#include <QVBoxLayout>
#include <QLabel>
#include <QFileDialog>
#include "cut.h"
#include "polyhedral_decomposition.h"
#include "subdivision_helper.h"
#include "export_helper.h"
#include "mesh_smoother.h"
#include "finalization.h"
#include "mesh_extractor.h"

extern bool batch_mode;

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void init_gui(GUI & gui, GlobalState & state)
{
    gui.but_cut_vol         = new QPushButton("Cut", &gui.window);
    gui.but_cut_conv        = new QPushButton("Cut To Conv", &gui.window);
    gui.but_polymesh        = new QPushButton("Meta Mesh", &gui.window);
    gui.but_finalize        = new QPushButton("Finalize MM", &gui.window);
    gui.but_subdivide       = new QPushButton("Subdivide", &gui.window);
    gui.but_export_meshes   = new QPushButton("Export", &gui.window);
    gui.cb_show_hide_loops  = new QCheckBox("Show Loops", &gui.window);
    gui.cb_show_active_loop = new QCheckBox("Show Active Loop", &gui.window);
    gui.cb_show_loop_normal = new QCheckBox("Show Loop Normals", &gui.window);
    gui.cb_show_inner_cuts  = new QCheckBox("Show Cuts", &gui.window);
    gui.but_undo            = new QPushButton("undo", &gui.window);
    gui.but_redo            = new QPushButton("redo", &gui.window);
    gui.but_reset           = new QPushButton("reset", &gui.window);
    gui.sb_curr_loop.setMinimum(-1);
    gui.sb_curr_loop.setValue(-1);
    gui.sb_curr_loop.setMaximum(state.loops.size()-1);
    gui.left.setMinimumSize(QSize(500,500));
    gui.right.setMinimumSize(QSize(500,500));
    gui.left.push_obj(&state.m_srf);
    gui.left.push_obj(&state.m_vol);
    gui.left.push_obj(&state.loops);
    gui.cb_show_hide_loops->setChecked(true);
    gui.cb_show_active_loop->setChecked(true);
    gui.cb_show_loop_normal->setChecked(false);
#ifndef TEST_SRF_ONLY
    state.m_srf.show_mesh(false);
#endif
    state.m_srf.show_wireframe(false);
    state.m_srf.show_mesh_flat();
    state.m_srf.show_marked_edge(false);
    state.m_vol.show_out_wireframe(false);
    state.m_vol.show_in_wireframe(false);
    state.m_vol.show_marked_face(false);
    state.m_vol.show_marked_edge(false);
    QWidget *toolbar = new QWidget(&gui.window);
    QVBoxLayout *l_toolbar = new QVBoxLayout();
    l_toolbar->addWidget(new QLabel("Loop (-1=none):"));
    l_toolbar->addWidget(&gui.sb_curr_loop);
    l_toolbar->addWidget(gui.but_cut_vol);
    l_toolbar->addWidget(gui.but_cut_conv);
    l_toolbar->addWidget(gui.but_polymesh);
    l_toolbar->addWidget(gui.but_finalize);
    l_toolbar->addWidget(gui.but_subdivide);
    l_toolbar->addWidget(gui.but_export_meshes);
    l_toolbar->addStretch();
    l_toolbar->addWidget(gui.cb_show_hide_loops);
    l_toolbar->addWidget(gui.cb_show_active_loop);
    l_toolbar->addWidget(gui.cb_show_loop_normal);
    l_toolbar->addWidget(gui.cb_show_inner_cuts);
    l_toolbar->addStretch();
    l_toolbar->addWidget(gui.but_undo);
    l_toolbar->addWidget(gui.but_redo);
    l_toolbar->addWidget(gui.but_reset);
    l_toolbar->addStretch();
    toolbar->setLayout(l_toolbar);
    toolbar->setMaximumWidth(150);
    QGridLayout *l_window = new QGridLayout();
    l_window->addWidget(&gui.left,0,0);
    l_window->addWidget(toolbar,0,1);
#ifndef DONT_SHOW_RIGHT_GUI
    l_window->addWidget(&gui.right,0,2);
#endif
    gui.window.setLayout(l_window);
    gui.window.show();

    // TOGGLE LOOP MATCHING VISUAL DEBUG
    gui.left.callback_key_press = [&](GLcanvas *c, QKeyEvent *e)
    {
        if(e->key()==Qt::Key_Space)
        {
            int lid = gui.sb_curr_loop.value();
            gui.sb_curr_loop.setValue(lid+1);
            gui.left.updateGL();
        }
        else if(e->key()==Qt::Key_Backspace)
        {
            int lid = gui.sb_curr_loop.value();
            gui.sb_curr_loop.setValue(lid-1);
            gui.left.updateGL();
        }
        else if(e->key()==Qt::Key_C)
        {
            emit(gui.but_cut_vol->clicked());
        }
        else if(e->key()==Qt::Key_S)
        {
            smoother(state.m_poly, state.m_srf);
            c->updateGL();
        }
    };

    // MOUSE PICKING (LEFT GUI)
    gui.left.callback_mouse_press = [&](GLcanvas *c, QMouseEvent *e)
    {
        // CMD + CLICK ON LOOP: do volumetric cut
        if(e->modifiers() == Qt::ControlModifier)
        {
            vec3d p;
            vec2i click(e->x(), e->y());
            if(!c->unproject(click, p)) return;
            uint lid = state.loops.pick_loop(p);
            state.loops.set_active_loop(lid, true);
            gui.sb_curr_loop.setValue(lid);
            emit(gui.but_cut_vol->clicked());
        }
        // SHIFT + CLICK ON LOOP: select loop
        else if(e->modifiers() == Qt::ShiftModifier)
        {
            vec3d p;
            vec2i click(e->x(), e->y());
            if(!c->unproject(click, p)) return;
            uint lid = state.loops.pick_loop(p);            
            gui.sb_curr_loop.setValue(lid);
            state.loops.set_active_loop(lid, true);
        }
    };

    // MOUSE PICKING (RIGHT GUI)
    gui.right.callback_mouse_press = [&](GLcanvas *c, QMouseEvent *e)
    {
        // SHIFT + CLICK ON LOOP: select/deselect polyhedron
        if(e->modifiers() == Qt::ShiftModifier)
        {
            vec3d p;
            vec2i click(e->x(), e->y());
            if(!c->unproject(click, p)) return;
            uint fid = state.m_poly.pick_face(p);
            uint pid = state.m_poly.adj_f2p(fid).front();
            if(state.m_poly.poly_data(pid).flags[HIDDEN]) pid = state.m_poly.adj_f2p(fid).back();
            static bool hide_others = false;
            hide_others = !hide_others;
            for(uint id=0; id<state.m_poly.num_polys(); ++id)
            {
                if(id!=pid) state.m_poly.poly_data(id).flags[HIDDEN] = !hide_others;
            }
            state.m_poly.updateGL();

            c->pop_all_markers();
            if(hide_others)
            {
                for(uint vid : state.m_poly.adj_p2v(pid))
                {
                    uint val = state.m_poly.poly_vert_valence(pid,vid);
                    if(val!=3) c->push_marker(state.m_poly.vert(vid), "val "+std::to_string(val), Color::BLUE());
                }
            }
            gui.right.updateGL();
        }
    };
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void init_events(GUI & gui, GlobalState & state)
{
    QSpinBox::connect(&gui.sb_curr_loop, static_cast<void(QSpinBox::*)(int)>(&QSpinBox::valueChanged), [&]()
    {
        int lid = gui.sb_curr_loop.value();
        switch(lid)
        {
            case -1 : state.loops.deactivate_all();          break;
            default : state.loops.set_active_loop(lid,true); break;
        }
        gui.left.updateGL();
    });

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    QPushButton::connect(gui.but_cut_vol, &QPushButton::clicked, [&]()
    {
        int lid = gui.sb_curr_loop.value();
        if(lid>=0)
        {
            cut(state, lid);
            emit(gui.but_polymesh->clicked());
            gui.sb_curr_loop.setValue(lid+1);
            state.m_vol.updateGL();
            gui.left.updateGL();
            if(!batch_mode) state.history.push_state(State(state.m_vol, state.loops));
        }
    });

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    QPushButton::connect(gui.but_cut_conv, &QPushButton::clicked, [&]()
    {
        bool converged = false;
        batch_mode = true;
        for(uint i=0; i<state.loops.size(); ++i)
        {
            if(state.loops.at(i).type==CONVEX) continue;

            if(state.loops.at(i).type==CONCAVE || !converged)
            {
                cut(state, i);
                state.profiler.push("Mesh Extractor");
                MeshExtractor me(state.m_vol);
                state.m_poly = me.mm;
                converged = me.converged(state.m_vol);
                state.profiler.pop("\n");
            }
        }
        batch_mode = false;
        state.m_vol.poly_color_wrt_label(false, 0.25, 1.0);
        state.m_vol.updateGL();
        gui.left.updateGL();

        emit(gui.but_polymesh->clicked());
    });

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    QPushButton::connect(gui.but_polymesh, &QPushButton::clicked, [&]()
    {        
        extract_polyhedral_mesh(state);
        gui.right.pop_all_occurrences_of(DRAWABLE_POLYHEDRALMESH);
        gui.right.push_obj(&state.m_poly);
        gui.right.updateGL();
    });

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    QPushButton::connect(gui.but_finalize, &QPushButton::clicked, [&]()
    {
        finalize_block_decomposition(state);
        classify_polyhedra(state.m_poly);
        gui.right.pop_all_occurrences_of(DRAWABLE_POLYHEDRALMESH);
        gui.right.push_obj(&state.m_poly);
        gui.right.updateGL();
    });

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    QPushButton::connect(gui.but_subdivide, &QPushButton::clicked, [&]()
    {
        SubdivisionHelper sh(state.m_vol, state.m_poly);
        state.m_poly = sh.subdivide();
        classify_polyhedra(state.m_poly);
        state.m_poly.updateGL();
        gui.right.updateGL();
    });

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    QPushButton::connect(gui.but_export_meshes, &QPushButton::clicked, [&]()
    {
        std::string filename = QFileDialog::getSaveFileName(NULL, "Export as Hexmesh", ".", "").toStdString();
        exporter(state, filename);
    });

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    QPushButton::connect(gui.cb_show_hide_loops, &QCheckBox::clicked, [&]()
    {
        state.loops.show_loops = gui.cb_show_hide_loops->isChecked();
        gui.left.updateGL();
    });

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    QPushButton::connect(gui.cb_show_active_loop, &QCheckBox::clicked, [&]()
    {
        state.loops.show_active_loop = gui.cb_show_active_loop->isChecked();
        gui.left.updateGL();
    });

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    QPushButton::connect(gui.cb_show_loop_normal, &QCheckBox::clicked, [&]()
    {
        state.loops.show_active_loop_normals = gui.cb_show_loop_normal->isChecked();
        gui.left.updateGL();
    });

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    QPushButton::connect(gui.cb_show_inner_cuts, &QCheckBox::clicked, [&]()
    {
        if(gui.cb_show_inner_cuts->isChecked())
        {
            state.m_vol.show_marked_edge(true);
            state.m_vol.show_marked_face(true);
            SlicerState ss;
            ss.X_thresh = 0.0;
            state.m_vol.slice(ss);
        }
        else
        {
            state.m_vol.show_marked_edge(false);
            state.m_vol.show_marked_face(false);
            SlicerState ss;
            state.m_vol.slice(ss);
        }
        gui.left.updateGL();
    });

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    QPushButton::connect(gui.but_undo, &QPushButton::clicked, [&]()
    {
        State tmp_state(state.m_vol, state.loops);
        if(state.history.undo(tmp_state))
        {
            state.m_vol = tmp_state.m_vol;
            state.loops = tmp_state.loops;
            state.m_srf.show_mesh(false);
            state.m_vol.show_mesh(true);
            state.m_poly.show_mesh(true);
            gui.left.updateGL();
        }
    });

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    QPushButton::connect(gui.but_redo, &QPushButton::clicked, [&]()
    {
        State tmp_state(state.m_vol, state.loops);
        if(state.history.redo(tmp_state))
        {
            state.m_vol = tmp_state.m_vol;
            state.loops = tmp_state.loops;
            state.m_srf.show_mesh(false);
            state.m_vol.show_mesh(true);
            state.m_poly.show_mesh(true);
            gui.left.updateGL();
        }
    });

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    QPushButton::connect(gui.but_reset, &QPushButton::clicked, [&]()
    {
        State obj(state.m_vol, state.loops);
        state.history.reset(obj);
        state.m_vol = obj.m_vol;
        state.loops = obj.loops;
        state.m_srf.show_mesh(false);
        state.m_vol.show_mesh(true);
        gui.left.updateGL();
    });
}
