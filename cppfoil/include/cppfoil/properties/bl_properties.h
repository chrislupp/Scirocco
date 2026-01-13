// cppfoil: C++ airfoil analysis library
// Copyright (c) 2024
// SPDX-License-Identifier: MIT

#ifndef CPPFOIL_PROPERTIES_BL_PROPERTIES_H
#define CPPFOIL_PROPERTIES_BL_PROPERTIES_H

#include "cppfoil/config.h"
#include "cppfoil/data/param.h"
#include "cppfoil/util/linearization.h"

namespace cppfoil {

// Shape parameter H = delta*/theta
ScalarLin4 GetH(const Vector4d& U);

// Kinematic shape parameter Hk (with compressibility correction)
ScalarLin4 GetHk(const Vector4d& U, const Param& param);

// Wake shape parameter with gap correction
ScalarLin4 GetHw(const Vector4d& U, double wgap);

// Kinetic energy shape parameter H*
ScalarLin4 GetHs(const Vector4d& U, const Param& param);

// Density shape parameter H**
ScalarLin4 GetHss(const Vector4d& U, const Param& param);

// Reynolds number based on momentum thickness
ScalarLin4 GetRet(const Vector4d& U, const Param& param);

// Skin friction coefficient
ScalarLin4 GetCf(const Vector4d& U, const Param& param);

// Dissipation coefficient
ScalarLin4 GetCDi(const Vector4d& U, const Param& param);

// Pressure coefficient from edge velocity
ScalarLin4 GetCpFromUe(double ue, const Param& param);

// Compressible velocity correction
ScalarLin4 GetUk(double ue, const Param& param);

// Squared Mach number
ScalarLin4 GetMach2(const Vector4d& U, const Param& param);

// Density from state
ScalarLin4 GetRho(const Vector4d& U, const Param& param);

// Amplification factor rate
ScalarLin4 GetDamp(const Vector4d& U, const Param& param);

// Initial ctau at transition
ScalarLin4 GetCttr(const Vector4d& U, const Param& param);

// Equilibrium ctau
ScalarLin4 GetCteq(const Vector4d& U, const Param& param);

// Upwind factor for stability
ScalarLin8 GetUpw(const Vector4d& U1, const Vector4d& U2, const Param& param);

// Upwind average of a quantity
ScalarLin8 Upwind(double upw, const Vector8d& upw_U,
                   double f1, const Vector4d& f1_U1,
                   double f2, const Vector4d& f2_U2);

}  // namespace cppfoil

#endif  // CPPFOIL_PROPERTIES_BL_PROPERTIES_H
