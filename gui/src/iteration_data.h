#pragma once

#include <QObject>
#include <QPointF>
#include <QString>
#include <QVector>

struct IterationData {
    Q_GADGET
    Q_PROPERTY(double alpha MEMBER alpha)
    Q_PROPERTY(double Re MEMBER Re)
    Q_PROPERTY(double Ma MEMBER Ma)
    Q_PROPERTY(double cl MEMBER cl)
    Q_PROPERTY(double cd MEMBER cd)
    Q_PROPERTY(double cdf MEMBER cdf)
    Q_PROPERTY(double cdp MEMBER cdp)
    Q_PROPERTY(double cdpi MEMBER cdpi)
    Q_PROPERTY(double cm MEMBER cm)
    Q_PROPERTY(bool converged MEMBER converged)
    Q_PROPERTY(double chord MEMBER chord)
    Q_PROPERTY(QString name MEMBER name)
    Q_PROPERTY(bool viscous MEMBER viscous)
    Q_PROPERTY(int foilN MEMBER foilN)

   public:
    double alpha{}, Re{}, Ma{};
    double cl{}, cd{}, cdf{}, cdp{}, cdpi{}, cm{};
    double chord{1.0};
    QString name;
    bool converged{}, viscous{};
    int foilN{};

    QVector<double> foil_x, foil_y;
    QVector<double> foil_tx, foil_ty;
    QVector<double> foil_s;
    QVector<double> wake_x, wake_y;
    int wakeN{};

    QVector<double> cp, cpi, cf, ue, uei, th, ds, sa, Ret, Hk;

    bool hasSensitivities{};
};
