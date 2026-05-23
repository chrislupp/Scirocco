// Compare cppfoil's R_V matrix entry-by-entry to mfoil's R_V at the
// same converged state (from comprehensive_reference.h5). Localizes which
// rows/columns of R_V differ, narrowing down where the 1e-5 lambda bug is.

#include <H5Cpp.h>
#include <Eigen/Dense>
#include <Eigen/SparseLU>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "cppfoil/core/adjoint.h"
#include "cppfoil/core/coupled_solver.h"
#include "cppfoil/core/forces.h"
#include "cppfoil/core/wake.h"
#include "cppfoil/cppfoil.h"

using namespace cppfoil;
using namespace Eigen;

static MatrixXd readMatrix(H5::H5File& f, const std::string& p) {
    H5::DataSet d = f.openDataSet(p);
    hsize_t dims[2];
    d.getSpace().getSimpleExtentDims(dims);
    std::vector<double> b(dims[0] * dims[1]);
    d.read(b.data(), H5::PredType::NATIVE_DOUBLE);
    MatrixXd m(dims[0], dims[1]);
    for (hsize_t i = 0; i < dims[0]; ++i)
        for (hsize_t j = 0; j < dims[1]; ++j) m(i, j) = b[i * dims[1] + j];
    return m;
}
static VectorXd readVector(H5::H5File& f, const std::string& p) {
    H5::DataSet d = f.openDataSet(p);
    hsize_t dims[1];
    d.getSpace().getSimpleExtentDims(dims);
    std::vector<double> b(dims[0]);
    d.read(b.data(), H5::PredType::NATIVE_DOUBLE);
    return Map<VectorXd>(b.data(), dims[0]);
}

int main() {
    const char* path = "/Users/clupp/tools/analysis/ProjectTrysil/cppfoil/"
                       "tests/comprehensive_reference.h5";
    H5::H5File f(path, H5F_ACC_RDONLY);
    MatrixXd foil_x = readMatrix(f, "comprehensive/foil_x");
    MatrixXd U_ref = readMatrix(f, "comprehensive/U");
    VectorXd wgap_ref = readVector(f, "comprehensive/wgap");
    MatrixXd R_V_ref;
    try {
        R_V_ref = readMatrix(f, "comprehensive/R_V");
    } catch (...) {
        std::cerr << "comprehensive/R_V not in HDF5; regenerate the reference\n";
        return 1;
    }

    int N = foil_x.cols();
    Cppfoil foil("0012", N - 1);
    foil.SetOper(3.0, -1.0, 1e6, 0.0);
    auto& s = foil.GetState();
    s.foil.x = foil_x;
    s.foil.N = N;
    RecomputeArcLength(s);
    s.oper.viscous = true;
    s.oper.initbl = true;
    s.param.verb = 0;
    s.isol.AIC.resize(0, 0);
    s.isol.gamref.resize(0, 0);
    s.isol.gam.resize(0);
    SolveViscous(s);
    s.glob.U = U_ref;
    s.vsol.wgap = wgap_ref;
    {
        VectorXd sstag_ue_ref = readVector(f, "comprehensive/sstag_ue");
        s.isol.sstag_ue = sstag_ue_ref;
    }
    s.glob.realloc = true;
    BuildGlobSys(s);
    StoreAdjointFactorization(s);

    int Nsys = s.glob.Nsys;
    int NN = 4 * Nsys;
    MatrixXd R_V_cpp = MatrixXd(s.glob.R_V);
    std::cout << "R_V cppfoil dims: " << R_V_cpp.rows() << "x"
              << R_V_cpp.cols() << "\n";
    std::cout << "R_V mfoil dims:   " << R_V_ref.rows() << "x"
              << R_V_ref.cols() << "\n";

    if (R_V_cpp.rows() != R_V_ref.rows()) return 1;

    // Per-row max abs error
    std::cout << "\nPer-row max abs error (top 20 rows by error):\n";
    std::vector<std::pair<double, int>> row_errs;
    for (int i = 0; i < NN; ++i) {
        double e = (R_V_cpp.row(i) - R_V_ref.row(i)).cwiseAbs().maxCoeff();
        row_errs.emplace_back(e, i);
    }
    std::sort(row_errs.begin(), row_errs.end(),
              [](auto a, auto b) { return a.first > b.first; });
    for (int k = 0; k < std::min<int>(20, NN); ++k) {
        int i = row_errs[k].second;
        const char* what = (i < 3 * Nsys) ? "R_BL"
                                            : "R_ue";
        int eq_idx = (i < 3 * Nsys) ? (i % 3) : -1;
        const char* eq[3] = {"mom", "shape", "lag"};
        int node = (i < 3 * Nsys) ? (i / 3) : (i - 3 * Nsys);
        std::cout << "  row " << i << " (" << what
                  << (eq_idx >= 0 ? std::string(" ") + eq[eq_idx] : "")
                  << " node=" << node << "): max abs="
                  << row_errs[k].first << "\n";
    }

    // Find worst single entry
    double max_abs = 0;
    int wi = 0, wj = 0;
    for (int i = 0; i < NN; ++i) {
        for (int j = 0; j < NN; ++j) {
            double e = std::abs(R_V_cpp(i, j) - R_V_ref(i, j));
            if (e > max_abs) {
                max_abs = e;
                wi = i;
                wj = j;
            }
        }
    }
    std::cout << "\nWorst entry: R_V[" << wi << "," << wj
              << "] cpp=" << R_V_cpp(wi, wj) << " mfoil=" << R_V_ref(wi, wj)
              << " abs=" << max_abs << "\n";

    // Stratify by column block (which variable of U[k,j] the column is)
    // col = 4*j + k where k in {0=th, 1=ds, 2=sa, 3=ue}
    std::cout << "\nMax abs error by column variable type:\n";
    const char* var[4] = {"th", "ds", "sa", "ue"};
    for (int k = 0; k < 4; ++k) {
        double mx = 0;
        int wj_in_block = 0;
        for (int j = 0; j < Nsys; ++j) {
            int col = 4 * j + k;
            double e =
                (R_V_cpp.col(col) - R_V_ref.col(col)).cwiseAbs().maxCoeff();
            if (e > mx) {
                mx = e;
                wj_in_block = j;
            }
        }
        std::cout << "  col " << var[k] << ": max=" << mx
                  << " at node " << wj_in_block << "\n";
    }

    // Count nonzero entries that differ significantly
    int n_diff = 0;
    int n_diff_bl = 0, n_diff_ue = 0;
    for (int i = 0; i < NN; ++i) {
        for (int j = 0; j < NN; ++j) {
            double e = std::abs(R_V_cpp(i, j) - R_V_ref(i, j));
            double mag =
                std::max(std::abs(R_V_cpp(i, j)), std::abs(R_V_ref(i, j)));
            if (mag > 1e-8 && e > 1e-10) {
                n_diff++;
                if (i < 3 * Nsys) n_diff_bl++;
                else n_diff_ue++;
            }
        }
    }
    std::cout << "\nEntries with |err|>1e-10 and |val|>1e-8: " << n_diff
              << " (BL rows: " << n_diff_bl
              << ", ue rows: " << n_diff_ue << ")\n";

    return 0;
}
