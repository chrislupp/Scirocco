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
#include <exception>
#include <iostream>
#include <string>
#include <vector>

#include "session.h"

#ifndef SCIROCCO_VERSION
#define SCIROCCO_VERSION "unknown"
#endif

/**
 * @brief fauXfoil: an XFOIL-style interactive front end to Scirocco.
 *
 * Reads commands from standard input, so it runs interactively or from a
 * script piped in ("fauXfoil < commands.txt"). An optional argument names an
 * airfoil file to load at startup.
 * @return Process exit status.
 */
int main(int argc, char** argv) {
    std::ios::sync_with_stdio(false);
    try {
        const std::vector<std::string> args(argv + 1, argv + argc);
        fauxfoil::Session session(std::cin, std::cout, SCIROCCO_VERSION);
        return session.Run(args);
    } catch (const std::exception& error) {
        std::cerr << "fauXfoil: " << error.what() << '\n';
        return 1;
    } catch (...) {
        return 1;
    }
}
