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
#include <QCursor>
#include <QGuiApplication>
#include <QImage>
#include <QKeyEvent>
#include <QMouseEvent>
#include <QPainter>
#include <QPen>
#include <QPointer>
#include <QScreen>
#include <QWidget>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <functional>
#include <future>
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

/// Where a scene's page lands in a target rectangle: the visible part of
/// the page (all of it, or the zoomed view) scaled to fit and centred.
struct PageMap {
    Rect view;
    double scale = 1.0;
    double ox = 0.0;    ///< Target x of view.x0
    double oy = 0.0;    ///< Target y of view.y0 (the bottom edge)
    double line = 1.0;  ///< Pen width

    PageMap(const Scene& scene, const QRectF& target) {
        view = scene.view.Empty()
                   ? Rect{0.0, 0.0, Scene::kWidth, Scene::kHeight}
                   : scene.view;
        const double w = view.x1 - view.x0;
        const double h = view.y1 - view.y0;
        scale = std::min(target.width() / w, target.height() / h);
        ox = target.left() + 0.5 * (target.width() - scale * w);
        oy = target.top() + 0.5 * (target.height() + scale * h);
        /* Lines stay thin when zoomed, as on XFOIL's screen */
        const double page_scale = std::min(target.width() / Scene::kWidth,
                                           target.height() / Scene::kHeight);
        line = std::max(1.0, 0.9 * page_scale);
    }
    [[nodiscard]] QPointF ToTarget(const Point& p) const {
        return {ox + scale * (p.x - view.x0), oy - scale * (p.y - view.y0)};
    }
    [[nodiscard]] Point ToPage(const QPointF& q) const {
        return {view.x0 + (q.x() - ox) / scale, view.y0 + (oy - q.y()) / scale};
    }
};

QPen PlotPen(const QColor& color, double width) {
    QPen pen(color, width);
    pen.setCapStyle(Qt::RoundCap);
    pen.setJoinStyle(Qt::RoundJoin);
    return pen;
}

}  // namespace

void RenderScene(QPainter& painter, const Scene& scene, const QRectF& target,
                 bool paper) {
    painter.save();
    painter.fillRect(target, paper ? Qt::white : Qt::black);
    painter.setRenderHint(QPainter::Antialiasing, true);
    painter.setClipRect(target);

    const PageMap map(scene, target);
    const double width = map.line;
    for (const Stroke& stroke : scene.strokes) {
        const QColor color = ColorOf(stroke.ink, paper);
        QPolygonF polygon;
        polygon.reserve(static_cast<int>(stroke.points.size()));
        for (const Point& p : stroke.points) polygon << map.ToTarget(p);
        if (!stroke.clip.Empty()) {
            const QPointF a = map.ToTarget({stroke.clip.x0, stroke.clip.y1});
            const QPointF b = map.ToTarget({stroke.clip.x1, stroke.clip.y0});
            painter.setClipRect(QRectF(a, b).intersected(target));
        }
        if (stroke.fill) {
            painter.setPen(Qt::NoPen);
            painter.setBrush(color);
            painter.drawPolygon(polygon);
        } else {
            QPen pen = PlotPen(color, width);
            const double unit = map.scale / width;
            if (stroke.dash == Dash::Dashed) {
                pen.setCapStyle(Qt::FlatCap);
                pen.setDashPattern({6.0 * unit, 4.0 * unit});
            } else if (stroke.dash == Dash::Dotted) {
                pen.setDashPattern({0.01, 4.0 * unit});
            }
            painter.setPen(pen);
            painter.setBrush(Qt::NoBrush);
            painter.drawPolyline(polygon);
        }
        if (!stroke.clip.Empty()) painter.setClipRect(target);
    }

    /* Lettering: stroke-font polylines with the same pen */
    for (const Label& label : scene.labels) {
        painter.setPen(PlotPen(ColorOf(label.ink, paper), width));
        painter.setBrush(Qt::NoBrush);
        for (const auto& stroke : LabelStrokes(label)) {
            QPolygonF polygon;
            polygon.reserve(static_cast<int>(stroke.size()));
            for (const Point& p : stroke) polygon << map.ToTarget(p);
            painter.drawPolyline(polygon);
        }
    }
    painter.restore();
}

namespace {

using Reply = std::function<void(std::optional<CursorHit>)>;

/// The plot window: draws the current scene scaled to fit, and takes
/// cursor input (clicks and keys) when the session asks for it.
class PlotWindow : public QWidget {
   public:
    PlotWindow() {
        setWindowTitle("fauXfoil");
        setAttribute(Qt::WA_ShowWithoutActivating);
        setAttribute(Qt::WA_OpaquePaintEvent);
        setFocusPolicy(Qt::StrongFocus);
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

    /// Wait for a click or key; @p reply is called once on the GUI thread.
    void BeginCursor(bool done_button, Reply reply) {
        Finish(std::nullopt);
        reply_ = std::move(reply);
        done_button_ = done_button;
        setCursor(Qt::CrossCursor);
        update();
    }

   protected:
    void paintEvent(QPaintEvent* /*event*/) override {
        QPainter painter(this);
        RenderScene(painter, scene_, QRectF(rect()), false);
        if (reply_ && done_button_) DrawDone(painter);
    }

    void mousePressEvent(QMouseEvent* event) override {
        if (!reply_) return;
        const QPointF at = event->position();
        CursorHit hit{PageMap(scene_, QRectF(rect())).ToPage(at), 0};
        if (done_button_ && DoneRect().contains(at)) hit.key = 'D';
        Finish(hit);
    }

    void keyPressEvent(QKeyEvent* event) override {
        if (!reply_ || event->text().isEmpty()) {
            QWidget::keyPressEvent(event);
            return;
        }
        const QPointF at = mapFromGlobal(QCursor::pos());
        const auto key = static_cast<char>(std::toupper(
            static_cast<unsigned char>(event->text()[0].toLatin1())));
        Finish(CursorHit{PageMap(scene_, QRectF(rect())).ToPage(at), key});
    }

    void closeEvent(QCloseEvent* event) override {
        Finish(std::nullopt);
        QWidget::closeEvent(event);
    }

   private:
    [[nodiscard]] QRectF DoneRect() const {
        const double w = 0.045 * width();
        const double h = 0.024 * height();
        return {width() - w - 6.0, height() - h - 6.0, w, h};
    }

    void DrawDone(QPainter& painter) const {
        const QRectF box = DoneRect();
        const QColor green(0, 200, 0);
        painter.setRenderHint(QPainter::Antialiasing, true);
        painter.setPen(PlotPen(green, 1.0));
        painter.setBrush(Qt::NoBrush);
        painter.drawRect(box);
        Label label;
        label.text = "Done";
        label.height = 0.5 * box.height();
        label.align = 0.5;
        label.at = {box.center().x(), -(box.center().y() + 0.5 * label.height)};
        for (const auto& stroke : LabelStrokes(label)) {
            QPolygonF polygon;
            for (const Point& p : stroke) polygon << QPointF(p.x, -p.y);
            painter.drawPolyline(polygon);
        }
    }

    void Finish(std::optional<CursorHit> hit) {
        if (!reply_) return;
        Reply reply = std::move(reply_);
        reply_ = nullptr;
        unsetCursor();
        update();
        reply(hit);
    }

    Scene scene_;
    Reply reply_;
    bool done_button_ = false;
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

std::optional<CursorHit> QtDisplay::Cursor(bool done_button) {
    auto promise = std::make_shared<std::promise<std::optional<CursorHit>>>();
    std::future<std::optional<CursorHit>> answer = promise->get_future();
    Impl* impl = impl_.get();
    QMetaObject::invokeMethod(
        qApp,
        [impl, promise, done_button]() noexcept {
            try {
                if (impl->window.isNull() || !impl->window->isVisible()) {
                    promise->set_value(std::nullopt);
                    return;
                }
                /* Cursor input wants the window's keys, so take focus
                   here (plot updates never do) */
                impl->window->raise();
                impl->window->activateWindow();
                impl->window->BeginCursor(
                    done_button, [promise](std::optional<CursorHit> hit) {
                        promise->set_value(hit);
                    });
            } catch (...) {
                std::fputs("fauXfoil: cursor input failed\n", stderr);
                try {
                    promise->set_value(std::nullopt);
                } catch (const std::future_error& error) {
                    /* Already answered by the window */
                    std::fputs(error.what(), stderr);
                }
            }
        },
        Qt::QueuedConnection);
    return answer.get();
}

}  // namespace fauxfoil
