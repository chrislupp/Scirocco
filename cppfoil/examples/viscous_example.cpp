// cppfoil example: Viscous NACA airfoil analysis
// Copyright (c) 2024
// SPDX-License-Identifier: MIT

#include <iostream>
#include <iomanip>
#include "cppfoil/cppfoil.h"

int main() {
    std::cout << "=== cppfoil Viscous Example: NACA 0012 at Re=1e6 ===" << std::endl;
    std::cout << std::endl;

    // Create a NACA 0012 airfoil with 129 panels
    cppfoil::Cppfoil foil("0012", 129);

    // Set Reynolds number for viscous analysis
    double Re = 1.0e6;
    double Ma = 0.0;  // Incompressible

    // Analyze at different angles of attack
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "  Alpha      CL        CD        CDf       CDp       CM" << std::endl;
    std::cout << "  -----   -------   -------   -------   -------   -------" << std::endl;

    for (double alpha = 0.0; alpha <= 10.0; alpha += 2.0) {
        // Set operating conditions (viscous)
        // SetOper(alpha, cl, Re, Ma) - use cl=-1 for given alpha
        foil.SetOper(alpha, -1.0, Re, Ma);

        // Solve
        bool converged = foil.Solve();

        // Get results
        double cl = foil.GetCl();
        double cd = foil.GetCd();
        double cdf = foil.GetCdf();
        double cdp = foil.GetCdp();
        double cm = foil.GetCm();

        // Print results
        std::cout << std::setw(7) << alpha << "  "
                  << std::setw(8) << cl << "  "
                  << std::setw(8) << cd << "  "
                  << std::setw(8) << cdf << "  "
                  << std::setw(8) << cdp << "  "
                  << std::setw(8) << cm;

        if (!converged) {
            std::cout << "  (not converged)";
        }
        std::cout << std::endl;

        // Debug: check if BL was computed
        if (alpha == 0.0) {
            const auto& state = foil.GetState();
            std::cout << "  Debug: Nsys=" << state.glob.Nsys
                      << ", U.cols=" << state.glob.U.cols()
                      << ", vsol.Is[0].size=" << state.vsol.Is[0].size()
                      << ", vsol.Is[1].size=" << state.vsol.Is[1].size()
                      << ", vsol.Is[2].size=" << state.vsol.Is[2].size()
                      << std::endl;
            if (state.glob.U.cols() > 0) {
                std::cout << "  Debug: th[0]=" << state.glob.U(0, 0)
                          << ", ds[0]=" << state.glob.U(1, 0)
                          << ", ue[0]=" << state.glob.U(3, 0)
                          << std::endl;
            }
        }
    }

    std::cout << std::endl;
    std::cout << "=== Analysis Complete ===" << std::endl;

    return 0;
}
