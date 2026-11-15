#include <scirocco/scirocco.h>

#include <cmath>

int main() {
    scirocco::Scirocco foil("0012", 60);
    foil.SetOper(4.0);
    return foil.Solve() && std::isfinite(foil.GetCl()) && foil.GetCl() > 0.4
               ? 0
               : 1;
}
