#include "hdf5_reader.h"

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace {

double readScalar(H5::Group& grp, const std::string& name) {
    auto ds = grp.openDataSet(name);
    double val{};
    ds.read(&val, H5::PredType::NATIVE_DOUBLE);
    return val;
}

int readInt(H5::Group& grp, const std::string& name) {
    auto ds = grp.openDataSet(name);
    int val{};
    ds.read(&val, H5::PredType::NATIVE_INT);
    return val;
}

QVector<double> readVector(H5::Group& grp, const std::string& name) {
    if (!grp.nameExists(name)) return {};
    auto ds = grp.openDataSet(name);
    auto space = ds.getSpace();
    hsize_t dims[1];
    space.getSimpleExtentDims(dims);
    QVector<double> vec(static_cast<int>(dims[0]));
    ds.read(vec.data(), H5::PredType::NATIVE_DOUBLE);
    return vec;
}

// Read [2, N] row-major matrix, return separate x and y vectors.
void readMatrix2xN(H5::Group& grp, const std::string& name,
                   QVector<double>& xout, QVector<double>& yout) {
    if (!grp.nameExists(name)) {
        xout.clear();
        yout.clear();
        return;
    }
    auto ds = grp.openDataSet(name);
    auto space = ds.getSpace();
    hsize_t dims[2];
    space.getSimpleExtentDims(dims);
    int rows = static_cast<int>(dims[0]);
    int cols = static_cast<int>(dims[1]);
    std::vector<double> buf(rows * cols);
    ds.read(buf.data(), H5::PredType::NATIVE_DOUBLE);
    xout.resize(cols);
    yout.resize(cols);
    for (int j = 0; j < cols; ++j) {
        xout[j] = buf[j];         // row 0
        yout[j] = buf[cols + j];  // row 1
    }
}

std::string readStringAttr(H5::Group& grp, const std::string& name) {
    if (!grp.attrExists(name)) return {};
    auto attr = grp.openAttribute(name);
    std::string val;
    H5::StrType str_type(H5::PredType::C_S1, H5T_VARIABLE);
    attr.read(str_type, val);
    return val;
}

}  // namespace

Hdf5Reader::Hdf5Reader(QObject* parent) : QObject(parent) {}

Hdf5Reader::~Hdf5Reader() = default;

void Hdf5Reader::openFile(const QString& path) {
    try {
        file_ =
            std::make_unique<H5::H5File>(path.toStdString(), H5F_ACC_RDONLY);
    } catch (H5::Exception&) {
        return;
    }

    auto attr = file_->openAttribute("num_iterations");
    attr.read(H5::PredType::NATIVE_INT, &numIterations_);

    filePath_ = path;
    currentIteration_ = 0;

    loadCoeffHistory();
    emit fileChanged();

    if (numIterations_ > 0) {
        loadIteration(0);
        buildSeries();
        emit iterationChanged();
    }
}

void Hdf5Reader::setCurrentIteration(int index) {
    if (index < 0 || index >= numIterations_ || index == currentIteration_)
        return;
    currentIteration_ = index;
    loadIteration(index);
    buildSeries();
    emit iterationChanged();
}

void Hdf5Reader::loadIteration(int index) {
    char grpName[32];
    std::snprintf(grpName, sizeof(grpName), "iteration_%03d", index);
    auto iter = file_->openGroup(grpName);

    auto oper = iter.openGroup("operating_conditions");
    currentData_.alpha = readScalar(oper, "alpha");
    currentData_.Re = readScalar(oper, "Re");
    currentData_.Ma = readScalar(oper, "Ma");
    currentData_.viscous = readInt(oper, "viscous") != 0;

    auto geom = iter.openGroup("geometry");
    currentData_.chord = readScalar(geom, "chord");
    currentData_.name = QString::fromStdString(readStringAttr(geom, "name"));

    auto panels = iter.openGroup("panels");
    auto foil = panels.openGroup("foil");
    currentData_.foilN = readInt(foil, "N");
    readMatrix2xN(foil, "x", currentData_.foil_x, currentData_.foil_y);
    readMatrix2xN(foil, "t", currentData_.foil_tx, currentData_.foil_ty);
    currentData_.foil_s = readVector(foil, "s");

    auto wake = panels.openGroup("wake");
    currentData_.wakeN = readInt(wake, "N");
    readMatrix2xN(wake, "x", currentData_.wake_x, currentData_.wake_y);

    auto coeff = iter.openGroup("coefficients");
    currentData_.cl = readScalar(coeff, "cl");
    currentData_.cd = readScalar(coeff, "cd");
    currentData_.cdf = readScalar(coeff, "cdf");
    currentData_.cdp = readScalar(coeff, "cdp");
    currentData_.cdpi = readScalar(coeff, "cdpi");
    currentData_.cm = readScalar(coeff, "cm");
    currentData_.converged = readInt(coeff, "converged") != 0;

    auto dist = iter.openGroup("distributions");
    currentData_.cp = readVector(dist, "cp");
    currentData_.cpi = readVector(dist, "cpi");
    currentData_.cf = readVector(dist, "cf");
    currentData_.ue = readVector(dist, "ue");
    currentData_.uei = readVector(dist, "uei");
    currentData_.th = readVector(dist, "th");
    currentData_.ds = readVector(dist, "ds");
    currentData_.sa = readVector(dist, "sa");
    currentData_.Ret = readVector(dist, "Ret");
    currentData_.Hk = readVector(dist, "Hk");

    currentData_.hasSensitivities =
        H5Lexists(iter.getId(), "sensitivities", H5P_DEFAULT) > 0;
}

void Hdf5Reader::loadCoeffHistory() {
    clHistory_.clear();
    cdHistory_.clear();
    cmHistory_.clear();

    for (int i = 0; i < numIterations_; ++i) {
        char grpName[32];
        std::snprintf(grpName, sizeof(grpName), "iteration_%03d", i);
        auto iter = file_->openGroup(grpName);
        auto coeff = iter.openGroup("coefficients");
        double cl = readScalar(coeff, "cl");
        double cd = readScalar(coeff, "cd");
        double cm = readScalar(coeff, "cm");
        clHistory_.append(QPointF(i, cl));
        cdHistory_.append(QPointF(i, cd));
        cmHistory_.append(QPointF(i, cm));
    }
}

void Hdf5Reader::buildSeries() {
    airfoilSeries_ = vecToSeries(currentData_.foil_x, currentData_.foil_y);
    wakeSeries_ = vecToSeries(currentData_.wake_x, currentData_.wake_y);

    int N = currentData_.foilN;
    double chord = currentData_.chord;

    // Cp vs x/c (foil portion only)
    cpSeries_.clear();
    cpiSeries_.clear();
    for (int i = 0; i < N && i < currentData_.cp.size(); ++i) {
        double xc = currentData_.foil_x[i] / chord;
        cpSeries_.append(QPointF(xc, currentData_.cp[i]));
        if (i < currentData_.cpi.size())
            cpiSeries_.append(QPointF(xc, currentData_.cpi[i]));
    }

    buildCpArrows();
}

void Hdf5Reader::buildCpArrows() {
    arrowBases_.clear();
    arrowTips_.clear();
    arrowCpValues_.clear();
    int N = currentData_.foilN;
    if (N == 0 || currentData_.cp.size() < N) return;

    int stride = (N > 150) ? 3 : (N > 80) ? 2 : 1;
    double scale = 0.1 * currentData_.chord;

    for (int i = 0; i < N; i += stride) {
        double tx = currentData_.foil_tx[i];
        double ty = currentData_.foil_ty[i];
        double mag = std::sqrt(tx * tx + ty * ty);
        if (mag < 1e-12) continue;
        double nx = -ty / mag;
        double ny = tx / mag;

        double cp = currentData_.cp[i];
        // Always point outward: use |Cp| along the outward normal
        double absCp = std::abs(cp);
        double dx = absCp * nx * scale;
        double dy = absCp * ny * scale;

        double bx = currentData_.foil_x[i];
        double by = currentData_.foil_y[i];
        arrowBases_.append(QPointF(bx, by));
        arrowTips_.append(QPointF(bx + dx, by + dy));
        arrowCpValues_.append(cp);
    }
}

QList<QPointF> Hdf5Reader::vecToSeries(const QVector<double>& xv,
                                       const QVector<double>& yv) const {
    QList<QPointF> pts;
    int n = std::min(xv.size(), yv.size());
    pts.reserve(n);
    for (int i = 0; i < n; ++i) pts.append(QPointF(xv[i], yv[i]));
    return pts;
}

QList<QPointF> Hdf5Reader::distributionSeries(const QString& name) const {
    int N = currentData_.foilN;
    const QVector<double>* dist = nullptr;

    if (name == "cf")
        dist = &currentData_.cf;
    else if (name == "ue")
        dist = &currentData_.ue;
    else if (name == "ue_inv")
        dist = &currentData_.uei;
    else if (name == "θ" || name == "th")
        dist = &currentData_.th;
    else if (name == "δ*" || name == "ds")
        dist = &currentData_.ds;
    else if (name == "n/Cτ" || name == "sa")
        dist = &currentData_.sa;
    else if (name == "Reθ" || name == "Ret")
        dist = &currentData_.Ret;
    else if (name == "Hk")
        dist = &currentData_.Hk;

    QList<QPointF> pts;
    if (!dist || dist->isEmpty()) return pts;
    double chord = currentData_.chord;
    for (int i = 0; i < N && i < dist->size(); ++i)
        pts.append(QPointF(currentData_.foil_x[i] / chord, (*dist)[i]));
    return pts;
}
