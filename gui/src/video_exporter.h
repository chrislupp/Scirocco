#pragma once

#include <QObject>
#include <QProcess>
#include <QQuickItem>
#include <QString>
#include <QTimer>

class Hdf5Reader;

class VideoExporter : public QObject {
    Q_OBJECT
    Q_PROPERTY(bool exporting READ exporting NOTIFY exportingChanged)
    Q_PROPERTY(int progress READ progress NOTIFY progressChanged)

   public:
    explicit VideoExporter(Hdf5Reader* reader, QObject* parent = nullptr);

    bool exporting() const { return exporting_; }
    int progress() const { return progress_; }

    Q_INVOKABLE void exportVideo(QQuickItem* captureItem,
                                 const QString& outputPath, bool darkMode,
                                 int fps = 10, int scale = 2);
    Q_INVOKABLE void cancel();

   signals:
    void exportingChanged();
    void progressChanged();
    void finished(const QString& path);
    void error(const QString& message);

   private slots:
    void captureNextFrame();

   private:
    Hdf5Reader* reader_;
    QQuickItem* captureItem_{};
    QProcess* ffmpeg_{};
    QTimer timer_;
    int totalFrames_{};
    int currentFrame_{};
    int originalIteration_{};
    bool targetDarkMode_{};
    bool exporting_{};
    int progress_{};
    int fps_{10};
    int grabW_{};
    int grabH_{};
    QString outputPath_;
};
