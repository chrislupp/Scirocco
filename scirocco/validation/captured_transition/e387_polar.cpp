/* Eppler 387 polar for the transition-length calibration.
 *
 * Usage: e387_polar COORDS RE L [NPANEL] [NCRIT] > out.txt
 *   COORDS  lower TE -> LE -> upper TE, one header line (e387_calibrate.py
 *           writes it); L = 0 selects the cut cell, L > 0 the captured
 *           model with capt_L = L. Prints alpha cl cd converged for
 *           alpha = -3..8 deg in 0.5 deg steps, each from a cold start.
 */
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "scirocco/scirocco.h"
using namespace scirocco;
int main(int argc, char** argv) {
  std::ifstream in(argv[1]); std::string line; std::getline(in, line);
  std::vector<double> v; double a, b;
  while (std::getline(in, line)) { std::istringstream ss(line); if (ss >> a >> b) { v.push_back(a); v.push_back(b); } }
  int n = v.size() / 2; Eigen::MatrixXd X(n, 2); for (int i = 0; i < n; ++i) { X(i,0) = v[2*i]; X(i,1) = v[2*i+1]; }
  double Re = atof(argv[2]), L = atof(argv[3]); int np = argc > 4 ? atoi(argv[4]) : 200; double ncrit = argc > 5 ? atof(argv[5]) : 9.0;
  for (double al = -3.0; al <= 8.01; al += 0.5) {
    Scirocco f(Eigen::MatrixXd(X.transpose()), np); auto& p = f.GetState().param; p.verb = 0; p.ncrit = ncrit;
    if (L > 0) { p.transition = TransitionModel::Captured; p.capt_L = L; }
    f.SetOper(al, -1.0, Re, 0.0);
    bool ok = false; try { ok = f.Solve(); } catch (...) {}
    std::printf("%.2f %.6f %.7f %d\n", al, ok ? f.GetCl() : 0.0, ok ? f.GetCd() : 0.0, (int)ok);
  }
}
