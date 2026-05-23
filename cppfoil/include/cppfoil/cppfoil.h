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
#ifndef CPPFOIL_CPPFOIL_H
#define CPPFOIL_CPPFOIL_H

#include <string>

#include "cppfoil/config.h"
#include "cppfoil/data/sensitivity.h"
#include "cppfoil/data/state.h"

namespace cppfoil {

// Main cppfoil solver class
class Cppfoil {
   public:
    // Construct from NACA digits
    explicit Cppfoil(const std::string& naca = "0012", int npanel = 199);

    // Construct from coordinates
    explicit Cppfoil(const MatrixXd& coords, int npanel = 199);

    // Set operating conditions
    // alpha : Angle of attack [degrees]
    // cl : Target lift coefficient (if set, alpha is computed)
    // Re : Reynolds number (enables viscous analysis)
    // Ma : Mach number
    void SetOper(double alpha = 0.0, double cl = -1.0, double Re = 0.0,
                 double Ma = 0.0);

    // Solve the current operating point
    // Returns true if converged (always true for inviscid)
    bool Solve();

    // Get results
    [[nodiscard]] double GetCl() const { return state_.post.cl; }
    [[nodiscard]] double GetCd() const { return state_.post.cd; }
    [[nodiscard]] double GetCdf() const { return state_.post.cdf; }
    [[nodiscard]] double GetCdp() const { return state_.post.cdp; }
    [[nodiscard]] double GetCm() const { return state_.post.cm; }
    [[nodiscard]] double GetAlpha() const { return state_.oper.alpha; }
    [[nodiscard]] bool IsConverged() const { return state_.glob.conv; }

    // Get distributions
    [[nodiscard]] const VectorXd& GetCp() const { return state_.post.cp; }
    [[nodiscard]] const VectorXd& GetCf() const { return state_.post.cf; }
    [[nodiscard]] const VectorXd& GetTheta() const { return state_.post.th; }
    [[nodiscard]] const VectorXd& GetDeltaStar() const {
        return state_.post.ds;
    }

    // Get panel coordinates
    [[nodiscard]] const MatrixXd& GetFoilCoords() const {
        return state_.foil.x;
    }
    [[nodiscard]] const MatrixXd& GetWakeCoords() const {
        return state_.wake.x;
    }

    // Geometry modification functions
    // Deploy a flap at hinge location with given deflection angle
    // xzhinge : Hinge location [x, z]
    // eta : Flap deflection angle in degrees (positive = down)
    void GeomFlap(const Vector2d& xzhinge, double eta);

    // Add camber to airfoil from given camber line increment
    // xzcamb : Camber increment points [2 x Nc] (x, z coords)
    void GeomAddCamber(const MatrixXd& xzcamb);

    // Derotate airfoil to make chord line horizontal
    void GeomDerotate();

    // Compute sensitivities of cl, cd, cm w.r.t. panel node coordinates
    // Must be called after Solve(). Uses adjoint method.
    // Returns Sensitivity struct with dcl_dx, dcd_dx, dcm_dx [2 x N] matrices
    Sensitivity ComputeSensitivities();

    // Accessors for sensitivity results (call ComputeSensitivities first)
    [[nodiscard]] const MatrixXd& GetDclDx() const {
        return sensitivity_.dcl_dx;
    }
    [[nodiscard]] const MatrixXd& GetDcdDx() const {
        return sensitivity_.dcd_dx;
    }
    [[nodiscard]] const MatrixXd& GetDcmDx() const {
        return sensitivity_.dcm_dx;
    }

    // Scalar operating-condition sensitivities
    [[nodiscard]] double GetDclDalpha() const {
        return sensitivity_.dcl_dalpha;
    }
    [[nodiscard]] double GetDcdDalpha() const {
        return sensitivity_.dcd_dalpha;
    }
    [[nodiscard]] double GetDcmDalpha() const {
        return sensitivity_.dcm_dalpha;
    }

    [[nodiscard]] double GetDclDMa() const { return sensitivity_.dcl_dMa; }
    [[nodiscard]] double GetDcdDMa() const { return sensitivity_.dcd_dMa; }
    [[nodiscard]] double GetDcmDMa() const { return sensitivity_.dcm_dMa; }

    [[nodiscard]] double GetDclDRe() const { return sensitivity_.dcl_dRe; }
    [[nodiscard]] double GetDcdDRe() const { return sensitivity_.dcd_dRe; }
    [[nodiscard]] double GetDcmDRe() const { return sensitivity_.dcm_dRe; }

    // Access internal state (for advanced use)
    FoilState& GetState() { return state_; }
    [[nodiscard]] const FoilState& GetState() const { return state_; }

   private:
    FoilState state_;
    Sensitivity sensitivity_;
};

}  // namespace cppfoil

#endif  // CPPFOIL_CPPFOIL_H
