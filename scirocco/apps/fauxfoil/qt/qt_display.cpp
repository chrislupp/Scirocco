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
#include "qt_display.h"

#include <QApplication>
#include <QGuiApplication>
#include <QImage>
#include <QPainter>
#include <QPainterPath>
#include <QPen>
#include <QPointer>
#include <QScreen>
#include <QWidget>
#include <cmath>
#include <cstdio>
#include <memory>

#include "stroke_font.h"

namespace fauxfoil {

namespace {

QColor ColorOf(Ink ink, bool paper) {
    switch (ink) {
        case Ink::Fore:
            return paper ? QColor(0, 0, 0) : QColor(255, 255, 255);
        case Ink::Red:
            return {255, 0, 0};
        case Ink::Blue:
            return {30, 140, 255};
        case Ink::DarkRed:
            return {170, 0, 0};
        case Ink::DarkBlue:
            return {20, 93, 170};
        case Ink::Green:
            return {0, 200, 0};
        case Ink::Orange:
            return {255, 140, 0};
        case Ink::Magenta:
            return {230, 0, 230};
        case Ink::Cyan:
            return {0, 210, 210};
        case Ink::Yellow:
            return {230, 230, 0};
    }
    return {255, 255, 255};
}

}  // namespace

void RenderScene(QPainter& painter, const Scene& scene, const QRectF& target,
                 bool paper) {
    painter.save();
    painter.fillRect(target, paper ? Qt::white : Qt::black);
    painter.setRenderHint(QPainter::Antialiasing, true);
    painter.setRenderHint(QPainter::TextAntialiasing, true);

    const double scale = std::min(target.width() / Scene::kWidth,
                                  target.height() / Scene::kHeight);
    const double ox =
        target.left() + 0.5 * (target.width() - scale * Scene::kWidth);
    const double oy =
        target.top() + 0.5 * (target.height() + scale * Scene::kHeight);
    auto map = [&](const Point& p) {
        return QPointF(ox + scale * p.x, oy - scale * p.y);
    };

    /* Lines a little under a page unit wide, never thinner than a pixel */
    const double width = std::max(1.0, 0.9 * scale);
    for (const Stroke& stroke : scene.strokes) {
        const QColor color = ColorOf(stroke.ink, paper);
        QPolygonF polygon;
        polygon.reserve(static_cast<int>(stroke.points.size()));
        for (const Point& p : stroke.points) polygon << map(p);
        if (stroke.fill) {
            painter.setPen(Qt::NoPen);
            painter.setBrush(color);
            painter.drawPolygon(polygon);
            continue;
        }
        QPen pen(color, width);
        pen.setCapStyle(Qt::RoundCap);
        pen.setJoinStyle(Qt::RoundJoin);
        if (stroke.dash == Dash::Dashed) {
            pen.setCapStyle(Qt::FlatCap);
            pen.setDashPattern({6.0 * scale / width, 4.0 * scale / width});
        } else if (stroke.dash == Dash::Dotted) {
            pen.setDashPattern({0.01, 4.0 * scale / width});
        }
        painter.setPen(pen);
        painter.setBrush(Qt::NoBrush);
        painter.drawPolyline(polygon);
    }

    /* Lettering: stroke-font polylines with the same pen */
    for (const Label& label : scene.labels) {
        QPen pen(ColorOf(label.ink, paper), width);
        pen.setCapStyle(Qt::RoundCap);
        pen.setJoinStyle(Qt::RoundJoin);
        painter.setPen(pen);
        painter.setBrush(Qt::NoBrush);
        for (const auto& stroke : LabelStrokes(label)) {
            QPolygonF polygon;
            polygon.reserve(static_cast<int>(stroke.size()));
            for (const Point& p : stroke) polygon << map(p);
            painter.drawPolyline(polygon);
        }
    }
    painter.restore();
}

namespace {

/// The plot window: draws the current scene scaled to fit.
class PlotWindow : public QWidget {
   public:
    PlotWindow() {
        setWindowTitle("fauXfoil");
        setAttribute(Qt::WA_ShowWithoutActivating);
        setAttribute(Qt::WA_OpaquePaintEvent);
        setMinimumSize(320, 247);
        QSize size(1245, 962);
        if (const QScreen* screen = QGuiApplication::primaryScreen()) {
            const QRect area = screen->availableGeometry();
            const double fit = std::min(0.8 * area.width() / 1245.0,
                                        0.8 * area.height() / 962.0);
            size = QSize(static_cast<int>(1245.0 * fit),
                         static_cast<int>(962.0 * fit));
        }
        resize(size);
    }

    void SetScene(Scene scene) {
        scene_ = std::move(scene);
        update();
    }

   protected:
    void paintEvent(QPaintEvent* /*event*/) override {
        QPainter painter(this);
        RenderScene(painter, scene_, QRectF(rect()), false);
    }

   private:
    Scene scene_;
};

}  // namespace

struct QtDisplay::Impl {
    QPointer<PlotWindow> window;
    std::string snapshot_dir;
    int snapshots = 0;

    void Show(const Scene& scene) {
        if (!snapshot_dir.empty()) Save(scene);
        if (window.isNull()) window = new PlotWindow;
        window->SetScene(scene);
        if (!window->isVisible()) window->show();
    }

    void Save(const Scene& scene) {
        constexpr double kRatio = 2.0;
        QImage image(static_cast<int>(Scene::kWidth * kRatio),
                     static_cast<int>(Scene::kHeight * kRatio),
                     QImage::Format_RGB32);
        image.setDevicePixelRatio(kRatio);
        QPainter painter(&image);
        RenderScene(painter, scene, QRectF(0, 0, Scene::kWidth, Scene::kHeight),
                    false);
        painter.end();
        char name[32];
        std::snprintf(name, sizeof name, "/plot%03d.png", ++snapshots);
        image.save(QString::fromStdString(snapshot_dir + name));
    }
};

QtDisplay::QtDisplay(std::string snapshot_dir)
    : impl_(std::make_unique<Impl>()) {
    impl_->snapshot_dir = std::move(snapshot_dir);
}

QtDisplay::~QtDisplay() {
    if (!impl_->window.isNull()) delete impl_->window.data();
}

void QtDisplay::Show(const Scene& scene) {
    Impl* impl = impl_.get();
    auto copy = std::make_shared<const Scene>(scene);
    QMetaObject::invokeMethod(
        qApp,
        [impl, copy]() noexcept {
            try {
                impl->Show(*copy);
            } catch (...) {
                std::fputs("fauXfoil: the plot could not be drawn\n", stderr);
            }
        },
        Qt::QueuedConnection);
}

}  // namespace fauxfoil
