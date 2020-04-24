#ifndef UNDOREDO_H
#define UNDOREDO_H

#include <cinolib/undo_redo.h>
#include "definitions.h"
#include "loops.h"

using namespace cinolib;

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

struct State
{
    State(const TetMesh & m_vol, const Loops & loops)
        : m_vol(m_vol)
        , loops(loops){}
    TetMesh  m_vol;
    Loops    loops;
};
typedef struct State State;

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

class UndoRedo : public AbstractUndoRedo<State,State>
{
    public:

        UndoRedo() : AbstractUndoRedo<State,State>() {}
        UndoRedo(const State & state) : AbstractUndoRedo<State,State>(state){}

        //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        void set_state(const State & prev, State & curr) const
        {
            curr.m_vol  = prev.m_vol;
            curr.loops  = prev.loops;
        }
};

#endif // UNDOREDO_H
