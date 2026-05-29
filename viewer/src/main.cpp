#include <QApplication>
#include <QQmlApplicationEngine>
#include <QQmlContext>
#include <QQuickStyle>

#include "hdf5_reader.h"
#include "video_exporter.h"

int main(int argc, char* argv[]) {
    QApplication app(argc, argv);
    app.setApplicationName("Scirocco Viewer");
    app.setOrganizationName("Scirocco");
    QQuickStyle::setStyle("Fusion");

    Hdf5Reader reader;
    VideoExporter exporter(&reader);

    QQmlApplicationEngine engine;
    engine.rootContext()->setContextProperty("reader", &reader);
    engine.rootContext()->setContextProperty("videoExporter", &exporter);
    engine.loadFromModule("SciroccoViewer", "Main");

    if (engine.rootObjects().isEmpty()) return -1;

    if (argc > 1) reader.openFile(QString(argv[1]));

    return app.exec();
}
