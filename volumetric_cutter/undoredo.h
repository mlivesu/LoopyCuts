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
