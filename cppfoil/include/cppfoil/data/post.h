// cppfoil: C++ airfoil analysis library
// Copyright (c) 2024
// SPDX-License-Identifier: MIT

#ifndef CPPFOIL_DATA_POST_H
#define CPPFOIL_DATA_POST_H

#include "cppfoil/config.h"

namespace cppfoil {

// Post-processing outputs and distributions
struct Post {
    VectorXd cp;              // cp distribution
    VectorXd cpi;             // Inviscid cp distribution
    double cl = 0.0;          // Lift coefficient
    VectorXd cl_ue;           // Linearization of cl w.r.t. ue [N, airfoil only]
    double cl_alpha = 0.0;    // Linearization of cl w.r.t. alpha
    double cm = 0.0;          // Moment coefficient
    double cdpi = 0.0;        // Near-field pressure drag coefficient
    double cd = 0.0;          // Total drag coefficient
    Vector4d cd_U;            // Linearization of cd w.r.t. last wake state [4]
    double cdf = 0.0;         // Skin friction drag coefficient
    double cdp = 0.0;         // Pressure drag coefficient

    // Distributions
    VectorXd th;              // Momentum thickness distribution
    VectorXd ds;              // Displacement thickness distribution
    VectorXd sa;              // Amplification factor/shear lag coeff distribution
    VectorXd ue;              // Edge velocity (compressible) distribution
    VectorXd uei;             // Inviscid edge velocity (compressible) distribution
    VectorXd cf;              // Skin friction distribution
    VectorXd Ret;             // Re_theta distribution
    VectorXd Hk;              // Kinematic shape parameter distribution
};

}  // namespace cppfoil

#endif  // CPPFOIL_DATA_POST_H
