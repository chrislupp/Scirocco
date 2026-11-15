/*
    Scirocco: C++ Airfoil Analysis Library

    Example: HDF5 output of a viscous alpha sweep.
    Each alpha point is stored as an iteration in a single HDF5 file.
*/
#include <cstdio>

#include "scirocco/io/hdf5_writer.h"
#include "scirocco/scirocco.h"

int main() {
    scirocco::Scirocco foil("2412", 199);
    scirocco::Hdf5Writer writer("alpha_sweep.h5");

    double Re = 1e6;
    double Ma = 0.0;

    for (double alpha = -2.0; alpha <= 10.0; alpha += 1.0) {
        foil.SetOper(alpha, -1.0, Re, Ma);
        bool converged = foil.Solve();

        if (converged) {
            auto sens = foil.ComputeSensitivities();
            writer.WriteIteration(foil, &sens);
        } else {
            writer.WriteIteration(foil);
        }

        std::printf("alpha=%5.1f  cl=%7.4f  cd=%8.5f  cm=%7.4f  %s\n", alpha,
                    foil.GetCl(), foil.GetCd(), foil.GetCm(),
                    converged ? "OK" : "NOT CONVERGED");
    }

    std::printf("\nWrote %d iterations to alpha_sweep.h5\n",
                writer.NumIterations());
    return 0;
}
