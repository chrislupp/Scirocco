#include <QApplication>
#include <QCommandLineParser>
#include <QIcon>
#include <QQmlApplicationEngine>
#include <QQmlContext>
#include <QQuickStyle>
#include <QQuickWindow>
#include <QTimer>

#include "hdf5_reader.h"
#include "video_exporter.h"

int main(int argc, char* argv[]) {
    QApplication app(argc, argv);
    app.setApplicationName("Scirocco Viewer");
    app.setOrganizationName("Scirocco");
    app.setWindowIcon(QIcon(":/qt/qml/SciroccoViewer/resources/app-icon.svg"));
    QQuickStyle::setStyle("Fusion");

    QCommandLineParser parser;
    parser.setApplicationDescription("View Scirocco HDF5 results.");
    parser.addHelpOption();
    parser.addPositionalArgument("file", "HDF5 file to open.", "[file]");
    const QCommandLineOption themeOption(
        "theme", "Start in the light or dark theme.", "light|dark", "dark");
    const QCommandLineOption screenshotOption(
        "screenshot",
        "Save the window to an image once it has drawn, then quit.", "image");
    parser.addOption(themeOption);
    parser.addOption(screenshotOption);
    parser.process(app);

    Hdf5Reader reader;
    VideoExporter exporter(&reader);

    QQmlApplicationEngine engine;
    engine.rootContext()->setContextProperty("reader", &reader);
    engine.rootContext()->setContextProperty("videoExporter", &exporter);
    engine.loadFromModule("SciroccoViewer", "Main");

    if (engine.rootObjects().isEmpty()) return -1;

    if (auto* theme =
            engine.singletonInstance<QObject*>("SciroccoViewer", "Theme")) {
        theme->setProperty("darkMode",
                           parser.value(themeOption) != QLatin1String("light"));
    }

    const QStringList files = parser.positionalArguments();
    if (!files.isEmpty()) reader.openFile(files.first());

    if (parser.isSet(screenshotOption)) {
        auto* window =
            qobject_cast<QQuickWindow*>(engine.rootObjects().first());
        const QString path = parser.value(screenshotOption);
        // Give the window time to lay out and draw before grabbing it.
        QTimer::singleShot(1500, &app, [window, path]() {
            const bool saved = window && window->grabWindow().save(path);
            QCoreApplication::exit(saved ? 0 : 1);
        });
    }

    return app.exec();
}
