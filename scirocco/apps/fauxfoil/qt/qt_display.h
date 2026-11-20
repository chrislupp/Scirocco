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

/// @file qt_display.h
/// The plot window, drawn with Qt.

#include <memory>
#include <string>

#include "display.h"

class QPainter;
class QRectF;

namespace fauxfoil {

/// Draw @p scene into @p target, keeping its aspect ratio. On screen the
/// background is black and the foreground white; @p paper swaps them.
void RenderScene(QPainter& painter, const Scene& scene, const QRectF& target,
                 bool paper);

/// Shows plots in one window. Construct and destroy it on the GUI thread;
/// Show() may be called from any thread.
class QtDisplay : public Display {
   public:
    /// @param snapshot_dir If not empty, every plot shown is also saved there
    ///        as plotNNN.png.
    explicit QtDisplay(std::string snapshot_dir = "");
    ~QtDisplay() override;
    QtDisplay(const QtDisplay&) = delete;
    QtDisplay& operator=(const QtDisplay&) = delete;

    void Show(const Scene& scene) override;

   private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

}  // namespace fauxfoil
