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

#ifndef LOOPS_H
#define LOOPS_H

#include <vector>
#include <cinolib/drawable_object.h>
#include "definitions.h"

class Loops : public DrawableObject, public std::vector<Loop>
{
    public:

        Loops(const char *filename = "", SrfMesh *m_srf = NULL, TetMesh *m_vol = NULL);

        //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        void       draw(const float) const;
        vec3d      scene_center() const { return m_vol->bbox().center(); }
        float      scene_radius() const { return m_vol->bbox().diag();   }
        ObjectType object_type()  const { return DRAWABLE_SKELETON;      }

        //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        void set_active_loop(const uint curr, const bool deactivate_others);
        void activate_all();
        void deactivate_all();

        //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        void close_open_creases();

        //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        void trace_loop(const uint           lid,
                        std::vector<uint>  & vids,
                        std::vector<vec3d> & points,
                        std::vector<vec3d> & srf_normals,
                        std::vector<vec3d> & cut_normals) const;

        //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        uint pick_loop(const vec3d & p) const;

        //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        bool show_loops               = true;
        bool show_active_loop         = true;
        bool show_active_loop_normals = false;

        //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        float length;
        float thick;

    private:

        const SrfMesh *m_srf;
              TetMesh *m_vol;
};

#endif // LOOPS_H
