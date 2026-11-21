#pragma once

#include <H5Cpp.h>

#include <QList>
#include <QObject>
#include <QPointF>
#include <QString>
#include <memory>

#include "iteration_data.h"

class Hdf5Reader : public QObject {
    Q_OBJECT

    Q_PROPERTY(int numIterations READ numIterations NOTIFY fileChanged)
    Q_PROPERTY(int currentIteration READ currentIteration WRITE
                   setCurrentIteration NOTIFY iterationChanged)
    Q_PROPERTY(IterationData data READ data NOTIFY iterationChanged)
    Q_PROPERTY(QString filePath READ filePath NOTIFY fileChanged)

    Q_PROPERTY(
        QList<QPointF> airfoilSeries READ airfoilSeries NOTIFY iterationChanged)
    Q_PROPERTY(
        QList<QPointF> wakeSeries READ wakeSeries NOTIFY iterationChanged)
    Q_PROPERTY(QList<QPointF> cpSeries READ cpSeries NOTIFY iterationChanged)
    Q_PROPERTY(QList<QPointF> cpiSeries READ cpiSeries NOTIFY iterationChanged)
    Q_PROPERTY(
        QList<QPointF> arrowBases READ arrowBases NOTIFY iterationChanged)
    Q_PROPERTY(QList<QPointF> arrowTips READ arrowTips NOTIFY iterationChanged)
    Q_PROPERTY(
        QList<qreal> arrowCpValues READ arrowCpValues NOTIFY iterationChanged)
    Q_PROPERTY(bool hasArrows READ hasArrows NOTIFY iterationChanged)

    Q_PROPERTY(QList<QPointF> clHistory READ clHistory NOTIFY fileChanged)
    Q_PROPERTY(QList<QPointF> cdHistory READ cdHistory NOTIFY fileChanged)
    Q_PROPERTY(QList<QPointF> cmHistory READ cmHistory NOTIFY fileChanged)

   public:
    explicit Hdf5Reader(QObject* parent = nullptr);
    ~Hdf5Reader() override;

    int numIterations() const { return numIterations_; }
    int currentIteration() const { return currentIteration_; }
    IterationData data() const { return currentData_; }
    QString filePath() const { return filePath_; }

    QList<QPointF> airfoilSeries() const { return airfoilSeries_; }
    QList<QPointF> wakeSeries() const { return wakeSeries_; }
    QList<QPointF> cpSeries() const { return cpSeries_; }
    QList<QPointF> cpiSeries() const { return cpiSeries_; }
    QList<QPointF> arrowBases() const { return arrowBases_; }
    QList<QPointF> arrowTips() const { return arrowTips_; }
    QList<qreal> arrowCpValues() const { return arrowCpValues_; }
    bool hasArrows() const { return !arrowBases_.isEmpty(); }

    QList<QPointF> clHistory() const { return clHistory_; }
    QList<QPointF> cdHistory() const { return cdHistory_; }
    QList<QPointF> cmHistory() const { return cmHistory_; }

    Q_INVOKABLE void openFile(const QString& path);
    Q_INVOKABLE void setCurrentIteration(int index);
    Q_INVOKABLE QList<QPointF> distributionSeries(const QString& name) const;

   signals:
    void fileChanged();
    void iterationChanged();

   private:
    void loadIteration(int index);
    void loadCoeffHistory();
    void buildSeries();
    void buildCpArrows();

    QList<QPointF> vecToSeries(const QVector<double>& xv,
                               const QVector<double>& yv) const;

    std::unique_ptr<H5::H5File> file_;
    int numIterations_{};
    int currentIteration_{};
    QString filePath_;
    IterationData currentData_;

    QList<QPointF> airfoilSeries_;
    QList<QPointF> wakeSeries_;
    QList<QPointF> cpSeries_;
    QList<QPointF> cpiSeries_;
    QList<QPointF> arrowBases_;
    QList<QPointF> arrowTips_;
    QList<qreal> arrowCpValues_;

    QList<QPointF> clHistory_;
    QList<QPointF> cdHistory_;
    QList<QPointF> cmHistory_;
};
