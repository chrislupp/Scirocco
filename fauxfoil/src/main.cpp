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
#include <atomic>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <string>
#include <vector>

#include "session.h"

#ifdef FAUXFOIL_WITH_QT
#include <QApplication>
#include <thread>

#include "qt/qt_display.h"
#endif

#ifndef SCIROCCO_VERSION
#define SCIROCCO_VERSION "unknown"
#endif

namespace {

#ifdef FAUXFOIL_WITH_QT
/// Whether a window can be opened. macOS always has a window server for a
/// logged-in user; elsewhere an X11 or Wayland display must be named.
/// Whether plots were turned off with FAUXFOIL_NOGRAPHICS.
bool GraphicsTurnedOff() {
    const char* off = std::getenv("FAUXFOIL_NOGRAPHICS");
    return off != nullptr && *off != '\0' && std::string(off) != "0";
}

bool HaveDisplay() {
    if (GraphicsTurnedOff()) return false;
    if (std::getenv("QT_QPA_PLATFORM") != nullptr) return true;
#ifdef __APPLE__
    return true;
#else
    return std::getenv("DISPLAY") != nullptr ||
           std::getenv("WAYLAND_DISPLAY") != nullptr;
#endif
}

/// Run the session on a worker thread while Qt runs the plot window on the
/// main thread, as macOS requires.
int RunWithGraphics(int argc, char** argv) {
    QApplication app(argc, argv);
    QApplication::setApplicationName("fauXFOIL");
    QApplication::setQuitOnLastWindowClosed(false);
    const char* snapshots = std::getenv("FAUXFOIL_SNAPSHOT_DIR");
    fauxfoil::QtDisplay display(snapshots != nullptr ? snapshots : "");

    /* Qt has removed its own options from argv */
    const std::vector<std::string> args(argv + 1, argv + argc);
    std::atomic<int> status{0};
    std::atomic<bool> finished{false};
    std::thread worker([&] {
        try {
            fauxfoil::Session session(std::cin, std::cout, SCIROCCO_VERSION);
            session.SetDisplay(&display);
            status = session.Run(args);
        } catch (const std::exception& error) {
            std::cerr << "fauXFOIL: " << error.what() << '\n';
            status = 1;
        } catch (...) {
            status = 1;
        }
        finished = true;
        QMetaObject::invokeMethod(&app, &QCoreApplication::quit,
                                  Qt::QueuedConnection);
    });
    QApplication::exec();
    if (!finished) {
        /* Quit from the window while the session waits for input */
        std::cout << '\n' << std::flush;
        std::_Exit(status);
    }
    worker.join();
    return status;
}
#endif

/// Why this run has no plot window, for the user; empty when they turned
/// plots off themselves.
std::string TextOnlyNotice() {
#ifdef FAUXFOIL_WITH_QT
    if (GraphicsTurnedOff()) return "";
    return "\n Text only: no display found, so plots are not shown.\n"
           " HARD still writes plots to plot.ps.\n";
#else
    return "\n Text only: this fauXFOIL was built without its plot window\n"
           " (configure with -DFAUXFOIL_GRAPHICS=ON and Qt 6 for plots).\n"
           " HARD still writes plots to plot.ps.\n";
#endif
}

}  // namespace

/**
 * @brief fauXFOIL: an XFOIL-style interactive front end to Scirocco.
 *
 * Reads commands from standard input, so it runs interactively or from a
 * script piped in ("fauXFOIL < commands.txt"). An optional argument names an
 * airfoil file to load at startup. Plots open in a window when fauXFOIL is
 * built with Qt and a display is available.
 * @return Process exit status.
 */
int main(int argc, char** argv) {
    std::ios::sync_with_stdio(false);
    try {
#ifdef FAUXFOIL_WITH_QT
        if (HaveDisplay()) return RunWithGraphics(argc, argv);
#endif
        const std::vector<std::string> args(argv + 1, argv + argc);
        fauxfoil::Session session(std::cin, std::cout, SCIROCCO_VERSION);
        session.SetNotice(TextOnlyNotice());
        return session.Run(args);
    } catch (const std::exception& error) {
        std::cerr << "fauXFOIL: " << error.what() << '\n';
        return 1;
    } catch (...) {
        return 1;
    }
}
