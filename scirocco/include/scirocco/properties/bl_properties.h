/*
    Scirocco: C++ Airfoil Analysis Library

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
    number: AFRL-2026-0405.

    The views expressed are those of the authors and do not reflect the
    official guidance or position of the United States Government, the
    Department of Defense or of the United States Air Force.

    Statement from DoD: The Appearance of external hyperlinks does not
    constitute endorsement by the United States Department of Defense (DoD) of
    the linked websites, of the information, products, or services contained
    therein. The DoD does not exercise any editorial, security, or other
    control over the information you may find at these locations.
*/
#pragma once

#include "scirocco/config.h"
#include "scirocco/data/param.h"
#include "scirocco/util/linearization.h"

namespace scirocco {

/** @brief Shape parameter H = delta* /theta
 *  @param U State vector [theta, ds, sa, ue]
 *  @return ScalarLin4 value and linearization w.r.t. U
 */
ScalarLin4 GetH(const Vector4d& U);

/** @brief Kinematic shape parameter Hk (with compressibility correction)
 *  @param U     State vector [theta, ds, sa, ue]
 *  @param param Flow parameter struct
 *  @return ScalarLin4 value and linearization w.r.t. U
 */
ScalarLin4 GetHk(const Vector4d& U, const Param& param);

/** @brief Wake shape parameter with gap correction
 *  @param U    State vector [theta, ds, sa, ue]
 *  @param wgap Wake gap height
 *  @return ScalarLin4 value and linearization w.r.t. U
 */
ScalarLin4 GetHw(const Vector4d& U, double wgap);

/** @brief Kinetic energy shape parameter H*
 *  @param U     State vector [theta, ds, sa, ue]
 *  @param param Flow parameter struct
 *  @return ScalarLin4 value and linearization w.r.t. U
 */
ScalarLin4 GetHs(const Vector4d& U, const Param& param);

/** @brief Density shape parameter H**
 *  @param U     State vector [theta, ds, sa, ue]
 *  @param param Flow parameter struct
 *  @return ScalarLin4 value and linearization w.r.t. U
 */
ScalarLin4 GetHss(const Vector4d& U, const Param& param);

/** @brief Reynolds number based on momentum thickness
 *  @param U     State vector [theta, ds, sa, ue]
 *  @param param Flow parameter struct
 *  @return ScalarLin4 value and linearization w.r.t. U
 */
ScalarLin4 GetRet(const Vector4d& U, const Param& param);

/** @brief Skin friction coefficient
 *  @param U     State vector [theta, ds, sa, ue]
 *  @param param Flow parameter struct
 *  @return ScalarLin4 value and linearization w.r.t. U
 */
ScalarLin4 GetCf(const Vector4d& U, const Param& param);

/** @brief Dissipation coefficient (combined turbulent / laminar with floor switch)
 *  @param U     State vector [theta, ds, sa, ue]
 *  @param param Flow parameter struct
 *  @return ScalarLin4 value and linearization w.r.t. U
 */
ScalarLin4 GetCDi(const Vector4d& U, const Param& param);

/** @brief Sub-components of GetCDi (exposed for the analytic Re chain)
 *  @param U     State vector [theta, ds, sa, ue]
 *  @param param Flow parameter struct
 *  @return ScalarLin4 value and linearization w.r.t. U
 */
ScalarLin4 GetCDi_Lam(const Vector4d& U, const Param& param);

/** @brief Turbulent wall contribution to dissipation coefficient
 *  @param U     State vector [theta, ds, sa, ue]
 *  @param param Flow parameter struct
 *  @return ScalarLin4 value and linearization w.r.t. U
 */
ScalarLin4 GetCDi_TurbWall(const Vector4d& U, const Param& param);

/** @brief Outer layer contribution to dissipation coefficient (turbulent)
 *  @param U     State vector [theta, ds, sa, ue]
 *  @param param Flow parameter struct
 *  @return ScalarLin4 value and linearization w.r.t. U
 */
ScalarLin4 GetCDi_Outer(const Vector4d& U, const Param& param);

/** @brief Laminar stress contribution to dissipation coefficient (turbulent)
 *  @param U     State vector [theta, ds, sa, ue]
 *  @param param Flow parameter struct
 *  @return ScalarLin4 value and linearization w.r.t. U
 */
ScalarLin4 GetCDi_LamStress(const Vector4d& U, const Param& param);

/** @brief Laminar wake dissipation contribution
 *  @param U      State vector [theta, ds, sa, ue]
 *  @param paramin Flow parameter struct (turb flag overridden to false internally)
 *  @return ScalarLin4 value and linearization w.r.t. U
 */
ScalarLin4 GetCDi_LamWake(const Vector4d& U, const Param& paramin);

/** @brief Pressure coefficient from edge velocity
 *  @param ue    Edge velocity
 *  @param param Flow parameter struct
 *  @return ScalarLin4 value and linearization w.r.t. U
 */
ScalarLin4 GetCpFromUe(double ue, const Param& param);

/** @brief Compressible velocity correction
 *  @param ue    Edge velocity
 *  @param param Flow parameter struct
 *  @return ScalarLin4 value and linearization w.r.t. U
 */
ScalarLin4 GetUk(double ue, const Param& param);

/** @brief Squared Mach number
 *  @param U     State vector [theta, ds, sa, ue]
 *  @param param Flow parameter struct
 *  @return ScalarLin4 value and linearization w.r.t. U
 */
ScalarLin4 GetMach2(const Vector4d& U, const Param& param);

/** @brief Density from state
 *  @param U     State vector [theta, ds, sa, ue]
 *  @param param Flow parameter struct
 *  @return ScalarLin4 value and linearization w.r.t. U
 */
ScalarLin4 GetRho(const Vector4d& U, const Param& param);

/** @brief Amplification factor rate
 *  @param U     State vector [theta, ds, sa, ue]
 *  @param param Flow parameter struct
 *  @return ScalarLin4 value and linearization w.r.t. U
 */
ScalarLin4 GetDamp(const Vector4d& U, const Param& param);

/** @brief Initial ctau at transition
 *  @param U     State vector [theta, ds, sa, ue]
 *  @param param Flow parameter struct
 *  @return ScalarLin4 value and linearization w.r.t. U
 */
ScalarLin4 GetCttr(const Vector4d& U, const Param& param);

/** @brief Equilibrium ctau
 *  @param U     State vector [theta, ds, sa, ue]
 *  @param param Flow parameter struct
 *  @return ScalarLin4 value and linearization w.r.t. U
 */
ScalarLin4 GetCteq(const Vector4d& U, const Param& param);

/** @brief Upwind factor for stability
 *  @param U1    Upstream state vector [theta, ds, sa, ue]
 *  @param U2    Downstream state vector [theta, ds, sa, ue]
 *  @param param Flow parameter struct
 *  @return ScalarLin8 value and linearization w.r.t. [U1; U2]
 */
ScalarLin8 GetUpw(const Vector4d& U1, const Vector4d& U2, const Param& param);

/** @brief Upwind average of a quantity
 *  @param upw   Upwind factor
 *  @param upw_U Linearization of upwind factor w.r.t. [U1; U2]
 *  @param f1    Quantity at upstream station
 *  @param f1_U1 Linearization of f1 w.r.t. U1
 *  @param f2    Quantity at downstream station
 *  @param f2_U2 Linearization of f2 w.r.t. U2
 *  @return ScalarLin8 blended value and linearization w.r.t. [U1; U2]
 */
ScalarLin8 Upwind(double upw, const Vector8d& upw_U, double f1,
                  const Vector4d& f1_U1, double f2, const Vector4d& f2_U2);

/** @brief BL thickness measure delta (delta* + weighted theta)
 *  @param U     State vector [theta, ds, sa, ue]
 *  @param param Flow parameter struct
 *  @return ScalarLin4 value and linearization w.r.t. U
 */
ScalarLin4 GetDe(const Vector4d& U, const Param& param);

/** @brief Normalized wall slip velocity
 *  @param U     State vector [theta, ds, sa, ue]
 *  @param param Flow parameter struct
 *  @return ScalarLin4 value and linearization w.r.t. U
 */
ScalarLin4 GetUs(const Vector4d& U, const Param& param);

/** @brief Result type for functions with x-derivative */
struct ScalarLinX {
    double value;
    Vector4d lin_U;
    double lin_x;
};

/** @brief cf * x / theta
 *  @param U     State vector [theta, ds, sa, ue]
 *  @param x     Streamwise coordinate
 *  @param param Flow parameter struct
 *  @return ScalarLinX value, linearization w.r.t. U, and d/dx
 */
ScalarLinX GetCfxt(const Vector4d& U, double x, const Param& param);

/** @brief cDi * x / theta (dissipation * x / theta)
 *  @param U     State vector [theta, ds, sa, ue]
 *  @param x     Streamwise coordinate
 *  @param param Flow parameter struct
 *  @return ScalarLinX value, linearization w.r.t. U, and d/dx
 */
ScalarLinX GetCDixt(const Vector4d& U, double x, const Param& param);

/** @brief Equilibrium pressure gradient 1/ue * due/dx
 *  @param ds    Displacement thickness
 *  @param ds_U  Linearization of ds w.r.t. [U1; U2]
 *  @param cf    Skin friction coefficient
 *  @param cf_U  Linearization of cf w.r.t. [U1; U2]
 *  @param Hk    Kinematic shape parameter
 *  @param Hk_U  Linearization of Hk w.r.t. [U1; U2]
 *  @param Ret   Reynolds number based on momentum thickness
 *  @param Ret_U Linearization of Ret w.r.t. [U1; U2]
 *  @param param Flow parameter struct
 *  @return ScalarLin8 value and linearization w.r.t. [U1; U2]
 */
ScalarLin8 GetUq(double ds, const Vector8d& ds_U, double cf,
                 const Vector8d& cf_U, double Hk, const Vector8d& Hk_U,
                 double Ret, const Vector8d& Ret_U, const Param& param);

namespace detail {

/** @brief Derivatives of thermodynamic reference quantities w.r.t. Mach number.
 *  Used to propagate dMa through BL closure functions.
 */
struct ThermoMaDeriv {
    double KTb_dMa = 0;  /**< d(KTb)/dMa */
    double KTl_dMa = 0;  /**< d(KTl)/dMa */
    double H0_dMa = 0;   /**< d(H0)/dMa */
    double mu0_dMa = 0;  /**< d(mu0)/dMa */
    double rho0_dMa = 0; /**< d(rho0)/dMa */
};

/** @brief Compute derivatives of all thermodynamic reference quantities w.r.t. Ma.
 *  @param param  Flow parameter struct
 *  @param Vinf   Freestream velocity
 *  @param rhoinf Freestream density
 *  @return ThermoMaDeriv struct populated with all reference quantity derivatives
 */
ThermoMaDeriv ComputeThermoMaDeriv(const Param& param, double Vinf,
                                   double rhoinf);

/** @brief d(Uk)/dMa — Karman-Tsien corrected edge velocity.
 *  @param ue    Edge velocity
 *  @param param Flow parameter struct
 *  @param d     Precomputed thermodynamic Ma-derivatives
 *  @return scalar d(Uk)/dMa
 */
double GetUk_dMa(double ue, const Param& param, const ThermoMaDeriv& d);

/** @brief d(Mach2)/dMa — squared local Mach number.
 *  @param U     State vector [theta, ds, sa, ue]
 *  @param param Flow parameter struct
 *  @param d     Precomputed thermodynamic Ma-derivatives
 *  @return scalar d(Mach2)/dMa
 */
double GetMach2_dMa(const Vector4d& U, const Param& param,
                    const ThermoMaDeriv& d);

/** @brief d(Hk)/dMa — kinematic shape parameter.
 *  @param U     State vector [theta, ds, sa, ue]
 *  @param param Flow parameter struct
 *  @param d     Precomputed thermodynamic Ma-derivatives
 *  @return scalar d(Hk)/dMa
 */
double GetHk_dMa(const Vector4d& U, const Param& param, const ThermoMaDeriv& d);

/** @brief d(Hss)/dMa — density shape parameter.
 *  @param U     State vector [theta, ds, sa, ue]
 *  @param param Flow parameter struct
 *  @param d     Precomputed thermodynamic Ma-derivatives
 *  @return scalar d(Hss)/dMa
 */
double GetHss_dMa(const Vector4d& U, const Param& param,
                  const ThermoMaDeriv& d);

/** @brief d(Ret)/dMa — momentum-thickness Reynolds number.
 *  @param U     State vector [theta, ds, sa, ue]
 *  @param param Flow parameter struct
 *  @param d     Precomputed thermodynamic Ma-derivatives
 *  @return scalar d(Ret)/dMa
 */
double GetRet_dMa(const Vector4d& U, const Param& param,
                  const ThermoMaDeriv& d);

/** @brief d(Hs)/dMa — kinetic energy shape parameter.
 *  @param U     State vector [theta, ds, sa, ue]
 *  @param param Flow parameter struct
 *  @param d     Precomputed thermodynamic Ma-derivatives
 *  @return scalar d(Hs)/dMa
 */
double GetHs_dMa(const Vector4d& U, const Param& param, const ThermoMaDeriv& d);

/** @brief d(Us)/dMa — normalized wall slip velocity.
 *  @param U     State vector [theta, ds, sa, ue]
 *  @param param Flow parameter struct
 *  @param d     Precomputed thermodynamic Ma-derivatives
 *  @return scalar d(Us)/dMa
 */
double GetUs_dMa(const Vector4d& U, const Param& param, const ThermoMaDeriv& d);

/** @brief d(Cf)/dMa — skin friction coefficient.
 *  @param U     State vector [theta, ds, sa, ue]
 *  @param param Flow parameter struct
 *  @param d     Precomputed thermodynamic Ma-derivatives
 *  @return scalar d(Cf)/dMa
 */
double GetCf_dMa(const Vector4d& U, const Param& param, const ThermoMaDeriv& d);

/** @brief d(CDi)/dMa — dissipation coefficient.
 *  @param U     State vector [theta, ds, sa, ue]
 *  @param param Flow parameter struct
 *  @param d     Precomputed thermodynamic Ma-derivatives
 *  @return scalar d(CDi)/dMa
 */
double GetCDi_dMa(const Vector4d& U, const Param& param,
                  const ThermoMaDeriv& d);

/** @brief d(Cf*xt)/dMa — Mach sensitivity of the cf·x/theta product.
 *  @param U     State vector [theta, ds, sa, ue]
 *  @param x     Streamwise coordinate
 *  @param param Flow parameter struct
 *  @param d     Precomputed thermodynamic Ma-derivatives
 *  @return scalar d(Cf*xt)/dMa
 */
double GetCfxt_dMa(const Vector4d& U, double x, const Param& param,
                   const ThermoMaDeriv& d);

/** @brief d(CDi*xt)/dMa — Mach sensitivity of the cDi·x/theta product.
 *  @param U     State vector [theta, ds, sa, ue]
 *  @param x     Streamwise coordinate
 *  @param param Flow parameter struct
 *  @param d     Precomputed thermodynamic Ma-derivatives
 *  @return scalar d(CDi*xt)/dMa
 */
double GetCDixt_dMa(const Vector4d& U, double x, const Param& param,
                    const ThermoMaDeriv& d);

}  /* namespace detail */

}  /* namespace scirocco */
