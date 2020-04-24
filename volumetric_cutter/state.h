#ifndef STATE_H
#define STATE_H

#include <QWidget>
#include <QSpinBox>
#include <QCheckBox>
#include <QPushButton>
#include <cinolib/gui/qt/glcanvas.h>
#include <cinolib/profiler.h>
#include "definitions.h"
#include "loops.h"
#include "undoredo.h"

typedef struct
{
    QWidget       window;              // mainwindow
    GLcanvas      left;                // contains model and loops
    GLcanvas      right;               // contains meta mesh
    QSpinBox      sb_curr_loop;        // id of current loop
    QPushButton  *but_cut_vol;         // performs volumetric cut
    QPushButton  *but_cut_conv;        // performs volumetric cuts until convergence
    QPushButton  *but_polymesh;        // extract the polyhedral mesh
    QPushButton  *but_finalize;        // post convergence analysis and val2 fixing
    QPushButton  *but_subdivide;       // apply midpoint subdivision
    QPushButton  *but_smooth;          // smooth the polymesh
    QPushButton  *but_export_meshes;   // export hexahedral elements
    QCheckBox    *cb_show_hide_loops;  // shows/hides candidate cutting loops
    QCheckBox    *cb_show_active_loop;
    QCheckBox    *cb_show_loop_normal;
    QCheckBox    *cb_show_inner_cuts;  // show inner faces belonging to a cut
    QPushButton  *but_undo;            // rollback to previous cut
    QPushButton  *but_redo;            // re-apply rolled back cut
    QPushButton  *but_reset;           // rollback all cuts
}
GUI;

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

typedef struct
{
    Loops    loops;
    SrfMesh  m_srf;
    TetMesh  m_vol;
    MetaMesh m_poly;
    Profiler profiler;
    UndoRedo history;
}
GlobalState;

#endif // STATE_H
