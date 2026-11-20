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
/* Renders fauXFOIL plots with Qt off screen (QT_QPA_PLATFORM=offscreen). */

#include <gtest/gtest.h>

#include <QApplication>
#include <QColor>
#include <QImage>
#include <QKeyEvent>
#include <QMouseEvent>
#include <QPainter>
#include <QWidget>
#include <chrono>
#include <future>
#include <thread>

#include "plot_scene.h"
#include "qt/qt_display.h"

namespace {

/// The application main() creates, for the tests that need one.
void App() { ASSERT_NE(QApplication::instance(), nullptr); }

QImage Render(const fauxfoil::Scene& scene, bool paper) {
    App();
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
    App();
    QImage image(2000, 962, QImage::Format_RGB32);
    QPainter painter(&image);
    fauxfoil::RenderScene(painter, scene, QRectF(0, 0, 2000, 962), false);
    painter.end();
    /* Page centred horizontally: nothing drawn left of its edge */
    const int left = static_cast<int>((2000 - fauxfoil::Scene::kWidth) / 2);
    EXPECT_EQ(image.pixelColor(left - 20, 961), QColor(Qt::black));
    EXPECT_NE(image.pixelColor(left + 20, 961), QColor(Qt::black));
}

/// Run the event loop until @p done or a timeout.
template <class Done>
bool Pump(Done done) {
    const auto end =
        std::chrono::steady_clock::now() + std::chrono::seconds(10);
    while (!done()) {
        QApplication::processEvents(QEventLoop::AllEvents, 20);
        if (std::chrono::steady_clock::now() > end) return false;
    }
    return true;
}

QWidget* PlotWindow() {
    for (QWidget* w : QApplication::topLevelWidgets())
        if (w->windowTitle() == "fauXFOIL" && w->isVisible()) return w;
    return nullptr;
}

TEST(FauxfoilQt, CursorReportsClicksAndKeysInPageCoordinates) {
    App();
    fauxfoil::QtDisplay display;
    fauxfoil::Scene scene;
    scene.Segment(0.0, 0.0, 100.0, 100.0);
    display.Show(scene);
    ASSERT_TRUE(Pump([] { return PlotWindow() != nullptr; }));
    QWidget* window = PlotWindow();
    window->resize(1245, 962);
    QApplication::processEvents();

    /* A click in the middle of the window is the middle of the page */
    auto click = std::async(std::launch::async,
                            [&display] { return display.Cursor(false); });
    ASSERT_TRUE(
        Pump([window] { return window->cursor().shape() == Qt::CrossCursor; }));
    const QPointF centre(window->width() / 2.0, window->height() / 2.0);
    QMouseEvent press(QEvent::MouseButtonPress, centre,
                      window->mapToGlobal(centre), Qt::LeftButton,
                      Qt::LeftButton, Qt::NoModifier);
    QApplication::sendEvent(window, &press);
    ASSERT_TRUE(Pump([&click] {
        return click.wait_for(std::chrono::seconds(0)) ==
               std::future_status::ready;
    }));
    const auto hit = click.get();
    ASSERT_TRUE(hit.has_value());
    EXPECT_EQ(hit->key, 0);
    EXPECT_NEAR(hit->page.x, 0.5 * fauxfoil::Scene::kWidth, 2.0);
    EXPECT_NEAR(hit->page.y, 0.5 * fauxfoil::Scene::kHeight, 2.0);

    /* A key press reports the key, upper case */
    auto key = std::async(std::launch::async,
                          [&display] { return display.Cursor(true); });
    ASSERT_TRUE(
        Pump([window] { return window->cursor().shape() == Qt::CrossCursor; }));
    QKeyEvent typed(QEvent::KeyPress, Qt::Key_D, Qt::NoModifier, "d");
    QApplication::sendEvent(window, &typed);
    ASSERT_TRUE(Pump([&key] {
        return key.wait_for(std::chrono::seconds(0)) ==
               std::future_status::ready;
    }));
    const auto pressed = key.get();
    ASSERT_TRUE(pressed.has_value());
    EXPECT_EQ(pressed->key, 'D');

    /* Closing the window ends cursor input */
    auto closed = std::async(std::launch::async,
                             [&display] { return display.Cursor(false); });
    ASSERT_TRUE(
        Pump([window] { return window->cursor().shape() == Qt::CrossCursor; }));
    window->close();
    ASSERT_TRUE(Pump([&closed] {
        return closed.wait_for(std::chrono::seconds(0)) ==
               std::future_status::ready;
    }));
    EXPECT_FALSE(closed.get().has_value());
}

}  // namespace

/* Qt needs its application object to outlive every widget and font, and
   to be destroyed before main() returns; gtest_main cannot arrange that.
   The tests draw off screen unless told otherwise, so they also run where
   there is no display, including when CMake lists them at build time. */
int main(int argc, char** argv) {
    if (qEnvironmentVariableIsEmpty("QT_QPA_PLATFORM"))
        qputenv("QT_QPA_PLATFORM", "offscreen");
    QApplication app(argc, argv);
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
