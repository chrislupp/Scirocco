#include "video_exporter.h"

#include <QFile>
#include <QImage>
#include <QQuickItemGrabResult>

#include "hdf5_reader.h"

VideoExporter::VideoExporter(Hdf5Reader* reader, QObject* parent)
    : QObject(parent), reader_(reader) {
    timer_.setSingleShot(true);
    connect(&timer_, &QTimer::timeout, this, &VideoExporter::captureNextFrame);
}

void VideoExporter::exportVideo(QQuickItem* captureItem,
                                const QString& outputPath, bool darkMode,
                                int fps, int scale) {
    if (exporting_ || !captureItem || reader_->numIterations() < 2) return;

    captureItem_ = captureItem;
    outputPath_ = outputPath;
    totalFrames_ = reader_->numIterations();
    currentFrame_ = 0;
    originalIteration_ = reader_->currentIteration();
    targetDarkMode_ = darkMode;
    fps_ = fps;

    int logicalW = static_cast<int>(captureItem->width());
    int logicalH = static_cast<int>(captureItem->height());
    grabW_ = ((logicalW * scale) / 2) * 2;
    grabH_ = ((logicalH * scale) / 2) * 2;

    exporting_ = true;
    progress_ = 0;
    emit exportingChanged();
    emit progressChanged();

    QString ffmpegBin = "ffmpeg";
    for (auto path : {"/opt/homebrew/bin/ffmpeg", "/usr/local/bin/ffmpeg"}) {
        if (QFile::exists(path)) {
            ffmpegBin = path;
            break;
        }
    }

    ffmpeg_ = new QProcess(this);
    ffmpeg_->setProcessChannelMode(QProcess::ForwardedErrorChannel);

    QStringList args;
    args << "-y"
         << "-f" << "rawvideo"
         << "-pix_fmt" << "rgba"
         << "-s" << QString("%1x%2").arg(grabW_).arg(grabH_) << "-r"
         << QString::number(fps) << "-i" << "pipe:0"
         << "-c:v" << "libx264"
         << "-pix_fmt" << "yuv420p"
         << "-preset" << "medium"
         << "-crf" << "15" << outputPath;

    ffmpeg_->start(ffmpegBin, args);
    if (!ffmpeg_->waitForStarted(3000)) {
        exporting_ = false;
        emit exportingChanged();
        emit error("Failed to start ffmpeg. Is it installed?");
        delete ffmpeg_;
        ffmpeg_ = nullptr;
        return;
    }

    reader_->setCurrentIteration(0);
    timer_.start(200);
}

void VideoExporter::cancel() {
    if (!exporting_) return;
    timer_.stop();
    if (ffmpeg_) {
        ffmpeg_->closeWriteChannel();
        ffmpeg_->waitForFinished(3000);
        ffmpeg_->deleteLater();
        ffmpeg_ = nullptr;
    }
    reader_->setCurrentIteration(originalIteration_);
    exporting_ = false;
    emit exportingChanged();
}

void VideoExporter::captureNextFrame() {
    if (!exporting_ || !captureItem_) return;

    QSize targetSize(grabW_, grabH_);
    auto result = captureItem_->grabToImage(targetSize);
    if (!result) {
        cancel();
        emit error("Failed to grab frame");
        return;
    }

    connect(
        result.data(), &QQuickItemGrabResult::ready, this, [this, result]() {
            QImage img = result->image();
            if (img.size() != QSize(grabW_, grabH_))
                img = img.scaled(grabW_, grabH_, Qt::IgnoreAspectRatio,
                                 Qt::SmoothTransformation);
            if (img.format() != QImage::Format_RGBA8888)
                img = img.convertToFormat(QImage::Format_RGBA8888);

            int rowBytes = grabW_ * 4;
            for (int y = 0; y < img.height(); y++)
                ffmpeg_->write(reinterpret_cast<const char*>(img.scanLine(y)),
                               rowBytes);

            currentFrame_++;
            progress_ = (currentFrame_ * 100) / totalFrames_;
            emit progressChanged();

            if (currentFrame_ >= totalFrames_) {
                ffmpeg_->closeWriteChannel();
                ffmpeg_->waitForFinished(10000);
                ffmpeg_->deleteLater();
                ffmpeg_ = nullptr;

                reader_->setCurrentIteration(originalIteration_);
                exporting_ = false;
                emit exportingChanged();
                emit finished(outputPath_);
            } else {
                reader_->setCurrentIteration(currentFrame_);
                timer_.start(200);
            }
        });
}
