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
/* Renders fauXfoil plots with Qt off screen (QT_QPA_PLATFORM=offscreen). */

#include <gtest/gtest.h>

#include <QColor>
#include <QImage>
#include <QPainter>

#include "plot_scene.h"
#include "qt/qt_display.h"

namespace {

QImage Render(const fauxfoil::Scene& scene, bool paper) {
    QImage image(static_cast<int>(fauxfoil::Scene::kWidth),
                 static_cast<int>(fauxfoil::Scene::kHeight),
                 QImage::Format_RGB32);
    QPainter painter(&image);
    fauxfoil::RenderScene(painter, scene,
                          QRectF(0, 0, image.width(), image.height()), paper);
    painter.end();
    return image;
}

/// Page point (y up) to image pixel (y down).
QColor At(const QImage& image, double x, double y) {
    return image.pixelColor(static_cast<int>(x),
                            static_cast<int>(fauxfoil::Scene::kHeight - y));
}

TEST(FauxfoilQt, DrawsLinesAndLabelsOnBlack) {
    fauxfoil::Scene scene;
    scene.Segment(100.0, 500.5, 900.0, 500.5, fauxfoil::Ink::Red);
    scene.Text(100.0, 100.0, 60.0, "IHI");
    const QImage screen = Render(scene, false);
    EXPECT_EQ(At(screen, 20.0, 20.0), QColor(Qt::black));
    EXPECT_EQ(At(screen, 500.0, 500.5), QColor(255, 0, 0));
    /* The stem of the first "I" runs up the middle of its cell */
    const double stem = 100.0 + 0.5 * fauxfoil::kAdvance * 60.0;
    EXPECT_EQ(At(screen, stem, 130.0), QColor(Qt::white));

    const QImage paper = Render(scene, true);
    EXPECT_EQ(At(paper, 20.0, 20.0), QColor(Qt::white));
    EXPECT_EQ(At(paper, stem, 130.0), QColor(Qt::black));
}

TEST(FauxfoilQt, KeepsAspectRatioInWideTargets) {
    fauxfoil::Scene scene;
    scene.Segment(0.0, 0.0, fauxfoil::Scene::kWidth, 0.0);
    QImage image(2000, 962, QImage::Format_RGB32);
    QPainter painter(&image);
    fauxfoil::RenderScene(painter, scene, QRectF(0, 0, 2000, 962), false);
    painter.end();
    /* Page centred horizontally: nothing drawn left of its edge */
    const int left = static_cast<int>((2000 - fauxfoil::Scene::kWidth) / 2);
    EXPECT_EQ(image.pixelColor(left - 20, 961), QColor(Qt::black));
    EXPECT_NE(image.pixelColor(left + 20, 961), QColor(Qt::black));
}

}  // namespace
