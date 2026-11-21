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
    const QCommandLineOption setOption(
        "set",
        "Set a property of the main window, such as activeSection=1 or "
        "showCpArrows=true. Repeat for more than one.",
        "name=value");
    parser.addOption(themeOption);
    parser.addOption(screenshotOption);
    parser.addOption(setOption);
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

    QObject* root = engine.rootObjects().first();
    for (const QString& assignment : parser.values(setOption)) {
        const qsizetype eq = assignment.indexOf('=');
        const QByteArray name = assignment.left(eq).toUtf8();
        if (eq <= 0 || root->property(name.constData()).isNull()) {
            qCritical("--set: no main window property in '%s'",
                      qPrintable(assignment));
            return 1;
        }
        root->setProperty(name.constData(), assignment.mid(eq + 1));
    }

    const QStringList files = parser.positionalArguments();
    if (!files.isEmpty()) reader.openFile(files.first());

    if (parser.isSet(screenshotOption)) {
        auto* window = qobject_cast<QQuickWindow*>(root);
        const QString path = parser.value(screenshotOption);
        // Give the window time to lay out and draw before grabbing it.
        QTimer::singleShot(1500, &app, [window, path]() {
            const bool saved = window && window->grabWindow().save(path);
            QCoreApplication::exit(saved ? 0 : 1);
        });
    }

    return app.exec();
}
