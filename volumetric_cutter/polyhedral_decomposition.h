#ifndef POLYHEDRAL_DECOMPOSITION_H
#define POLYHEDRAL_DECOMPOSITION_H

#include "state.h"
#include "definitions.h"

void extract_polyhedral_mesh(GlobalState & state);

void classify_polyhedra(MetaMesh & m);

#endif // POLYHEDRAL_DECOMPOSITION_H
