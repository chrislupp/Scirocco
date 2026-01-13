// cppfoil: C++ airfoil analysis library
// Copyright (c) 2024
// SPDX-License-Identifier: MIT

#ifndef CPPFOIL_CPPFOIL_H
#define CPPFOIL_CPPFOIL_H

#include <string>
#include "cppfoil/config.h"
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
    void SetOper(double alpha = 0.0, double cl = -1.0,
                 double Re = 0.0, double Ma = 0.0);

    // Solve the current operating point
    // Returns true if converged (always true for inviscid)
    bool Solve();

    // Get results
    double GetCl() const { return state_.post.cl; }
    double GetCd() const { return state_.post.cd; }
    double GetCdf() const { return state_.post.cdf; }
    double GetCdp() const { return state_.post.cdp; }
    double GetCm() const { return state_.post.cm; }
    double GetAlpha() const { return state_.oper.alpha; }
    bool IsConverged() const { return state_.glob.conv; }

    // Get distributions
    const VectorXd& GetCp() const { return state_.post.cp; }
    const VectorXd& GetCf() const { return state_.post.cf; }
    const VectorXd& GetTheta() const { return state_.post.th; }
    const VectorXd& GetDeltaStar() const { return state_.post.ds; }

    // Get panel coordinates
    const MatrixXd& GetFoilCoords() const { return state_.foil.x; }
    const MatrixXd& GetWakeCoords() const { return state_.wake.x; }

    // Access internal state (for advanced use)
    FoilState& GetState() { return state_; }
    const FoilState& GetState() const { return state_; }

private:
    FoilState state_;
};

}  // namespace cppfoil

#endif  // CPPFOIL_CPPFOIL_H
