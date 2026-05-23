// Compare cppfoil analytic ∂R_BL/∂Ma vs mfoil-CS reference
// (drv_dma_cs_reference.h5). CS-based — replaces FD-based probe_drv_dma.cpp.
#include <H5Cpp.h>

#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <vector>

#include "cppfoil/core/coupled_solver.h"
#include "cppfoil/core/inviscid.h"
#include "cppfoil/cppfoil.h"

using namespace cppfoil;
using namespace Eigen;

int main() {
    const char* h5path =
        "/Users/clupp/tools/analysis/ProjectTrysil/cppfoil/tests/"
        "drv_dma_cs_reference.h5";
    H5::H5File file(h5path, H5F_ACC_RDONLY);
    int Nsys = 0;
    file.openDataSet("Nsys").read(&Nsys, H5::PredType::NATIVE_INT);
    std::vector<double> dR_cs_vec(3 * Nsys);
    file.openDataSet("dR_dMa").read(dR_cs_vec.data(),
                                    H5::PredType::NATIVE_DOUBLE);
    VectorXd dR_cs = Map<VectorXd>(dR_cs_vec.data(), 3 * Nsys);

    Cppfoil foil("0012", 100);
    foil.SetOper(3.0, -1.0, 1e6, 0.3);
    auto& s = foil.GetState();
    s.oper.viscous = true;
    s.oper.initbl = true;
    s.param.verb = 0;
    SolveViscous(s);
    if (s.glob.Nsys != Nsys) {
        std::cerr << "Nsys mismatch: cppfoil=" << s.glob.Nsys << " h5=" << Nsys
                  << "\n";
        return 1;
    }
    BuildGlobSys(s);
    VectorXd dR_an = s.glob.dR_dMa_times_Ma / s.oper.Ma;

    double max_abs = 0;
    int wi = 0;
    for (int i = 0; i < 3 * Nsys; ++i) {
        double e = std::abs(dR_an(i) - dR_cs(i));
        if (e > max_abs) {
            max_abs = e;
            wi = i;
        }
    }
    std::cout << "Worst: row " << wi << " station " << wi / 3 << " eq "
              << wi % 3 << " analytic=" << dR_an(wi) << " CS=" << dR_cs(wi)
              << " abs=" << max_abs << "\n";

    int n_diff = 0;
    std::cout << "Entries with abs > 1e-9 (showing up to 20):\n";
    for (int i = 0; i < 3 * Nsys; ++i) {
        double e = std::abs(dR_an(i) - dR_cs(i));
        if (e > 1e-9 && std::abs(dR_cs(i)) > 1e-10) {
            if (n_diff < 20) {
                double rel = e / std::max(std::abs(dR_cs(i)), 1e-12);
                std::cout << "  station " << i / 3 << " eq " << i % 3
                          << ": an=" << dR_an(i) << " cs=" << dR_cs(i)
                          << " abs=" << e << " rel=" << rel << "\n";
            }
            ++n_diff;
        }
    }
    std::cout << "Total entries with abs>1e-9: " << n_diff << "\n";
    std::cout << "Total |dR_an - dR_cs| L2: " << (dR_an - dR_cs).norm() << "\n";
    return 0;
}
