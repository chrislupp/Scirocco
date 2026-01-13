// cppfoil example: NACA airfoil analysis
// Copyright (c) 2024
// SPDX-License-Identifier: MIT

#include <iostream>
#include <iomanip>
#include "cppfoil/cppfoil.h"

int main() {
    std::cout << "=== cppfoil Example: NACA Airfoil Analysis ===" << std::endl;
    std::cout << std::endl;

    // Create a NACA 0012 airfoil with 199 panels
    cppfoil::Cppfoil foil("0012", 199);

    // Analyze at different angles of attack
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "  Alpha      CL        CD        CM" << std::endl;
    std::cout << "  -----   -------   -------   -------" << std::endl;

    for (double alpha = -4.0; alpha <= 12.0; alpha += 2.0) {
        // Set operating conditions (inviscid for now)
        foil.SetOper(alpha, -1.0, 0.0, 0.0);

        // Solve
        foil.Solve();

        // Print results
        std::cout << std::setw(7) << alpha << "  "
                  << std::setw(8) << foil.GetCl() << "  "
                  << std::setw(8) << foil.GetCd() << "  "
                  << std::setw(8) << foil.GetCm() << std::endl;
    }

    std::cout << std::endl;
    std::cout << "=== Analysis Complete ===" << std::endl;

    return 0;
}
