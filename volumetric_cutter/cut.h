#ifndef CUT_H
#define CUT_H

#include "state.h"

void cut(GlobalState & state, const uint lid);

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void revert_cut(GlobalState & state, const uint lid);

#endif // CUT_H
