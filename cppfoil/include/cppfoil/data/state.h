// cppfoil: C++ airfoil analysis library
// Copyright (c) 2024
// SPDX-License-Identifier: MIT

#ifndef CPPFOIL_DATA_STATE_H
#define CPPFOIL_DATA_STATE_H

#include <string>
#include "cppfoil/data/geom.h"
#include "cppfoil/data/panel.h"
#include "cppfoil/data/oper.h"
#include "cppfoil/data/isol.h"
#include "cppfoil/data/vsol.h"
#include "cppfoil/data/glob.h"
#include "cppfoil/data/post.h"
#include "cppfoil/data/param.h"

namespace cppfoil {

// Main state container for cppfoil solver
struct FoilState {
    std::string version = "1.0.0";
    Geom geom;            // Geometry
    Panel foil;           // Airfoil panels
    Panel wake;           // Wake panels
    Oper oper;            // Operating conditions
    Isol isol;            // Inviscid solution variables
    Vsol vsol;            // Viscous solution variables
    Glob glob;            // Global system variables
    Post post;            // Post-processing quantities
    Param param;          // Parameters
};

}  // namespace cppfoil

#endif  // CPPFOIL_DATA_STATE_H
