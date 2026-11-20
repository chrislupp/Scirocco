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

/// @file display.h
/// Where a session's plots go.

#include <optional>

#include "plot_scene.h"

namespace fauxfoil {

/// A mouse click or key press in the plot window.
struct CursorHit {
    Point page;    ///< Where, in page coordinates
    char key = 0;  ///< Upper-case key typed, or 0 for a mouse button
};

/// A plot window. Its methods may be called from any thread.
class Display {
   public:
    virtual ~Display() = default;
    /// Replace the plot on screen with @p scene.
    virtual void Show(const Scene& scene) = 0;
    /// Wait for a mouse click or key press in the window, as XFOIL's cursor
    /// input does. With @p done_button a "Done" button is shown; clicking
    /// it reports the key 'D'.
    /// @return Nothing when there is no window to click in.
    virtual std::optional<CursorHit> Cursor(bool done_button) {
        (void)done_button;
        return std::nullopt;
    }
};

}  // namespace fauxfoil
