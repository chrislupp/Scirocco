/*
    cppfoil: C++ Airfoil Analysis Library

    Copyright 2024-2026 Christopher A. Lupp

    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

        http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software
    distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and
    limitations under the License.

    This work has been cleared for public release, distribution unlimited, case
    number: AFRL-2026-XXXX.

    The views expressed are those of the authors and do not reflect the
    official guidance or position of the United States Government, the
    Department of Defense or of the United States Air Force.

    Statement from DoD: The Appearance of external hyperlinks does not
    constitute endorsement by the United States Department of Defense (DoD) of
    the linked websites, of the information, products, or services contained
    therein. The DoD does not exercise any editorial, security, or other
    control over the information you may find at these locations.
*/
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
ScalarLin8 Upwind(double upw, const Vector8d& upw_U, double f1,
                  const Vector4d& f1_U1, double f2, const Vector4d& f2_U2);

// BL thickness measure delta (delta* + weighted theta)
ScalarLin4 GetDe(const Vector4d& U, const Param& param);

// Normalized wall slip velocity
ScalarLin4 GetUs(const Vector4d& U, const Param& param);

// Result type for functions with x-derivative
struct ScalarLinX {
    double value;
    Vector4d lin_U;
    double lin_x;
};

// cf * x / theta
ScalarLinX GetCfxt(const Vector4d& U, double x, const Param& param);

// cDi * x / theta (dissipation * x / theta)
ScalarLinX GetCDixt(const Vector4d& U, double x, const Param& param);

// Equilibrium pressure gradient 1/ue * due/dx
ScalarLin8 GetUq(double ds, const Vector8d& ds_U, double cf,
                 const Vector8d& cf_U, double Hk, const Vector8d& Hk_U,
                 double Ret, const Vector8d& Ret_U, const Param& param);

}  // namespace cppfoil

#endif  // CPPFOIL_PROPERTIES_BL_PROPERTIES_H
