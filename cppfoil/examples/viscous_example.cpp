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

#include <iomanip>
#include <iostream>

#include "cppfoil/cppfoil.h"

int main() {
    std::cout << "=== cppfoil Viscous Example: NACA 0012 at Re=1e6 ==="
              << std::endl;
    std::cout << std::endl;

    // Create a NACA 0012 airfoil with 129 panels
    cppfoil::Cppfoil foil("0012", 129);

    // Set Reynolds number for viscous analysis
    double Re = 1.0e6;
    double Ma = 0.0;  // Incompressible

    // Analyze at different angles of attack
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "  Alpha      CL        CD        CDf       CDp       CM"
              << std::endl;
    std::cout << "  -----   -------   -------   -------   -------   -------"
              << std::endl;

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
        std::cout << std::setw(7) << alpha << "  " << std::setw(8) << cl << "  "
                  << std::setw(8) << cd << "  " << std::setw(8) << cdf << "  "
                  << std::setw(8) << cdp << "  " << std::setw(8) << cm;

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
                          << ", ue[0]=" << state.glob.U(3, 0) << std::endl;
            }
        }
    }

    std::cout << std::endl;
    std::cout << "=== Analysis Complete ===" << std::endl;

    return 0;
}
