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
*/

// The program shown on the documentation quickstart page.

#include <iostream>

#include "scirocco/scirocco.h"

int main() {
    scirocco::Scirocco foil("2412", 160);  // NACA 2412, 160 panels

    // alpha = 4 deg, no lift target, Re = 1e6, Mach 0
    foil.SetOper(4.0, -1.0, 1.0e6, 0.0);
    if (!foil.Solve()) {
        std::cerr << "did not converge\n";
        return 1;
    }

    std::cout << "Cl = " << foil.GetCl() << '\n'
              << "Cd = " << foil.GetCd() << '\n'
              << "Cm = " << foil.GetCm() << '\n';
    return 0;
}
