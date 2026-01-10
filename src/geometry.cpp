// Copyright (c) 2024 ProjectTrysil Authors
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Geometry utilities - explicit template instantiations.

#include "xfoil/geometry.h"

namespace xfoil {

// Explicit template instantiations for double
template NacaResult<double> naca4<double>(int, int);
template NacaResult<double> naca4_params<double>(double, double, double, int);
template NacaResult<double> naca5<double>(int, int);
template double lefind<double>(const VectorX<double>&, const VectorX<double>&,
                               const VectorX<double>&, const VectorX<double>&,
                               const VectorX<double>&);
template GeomParams<double> geopar<double>(const VectorX<double>&,
                                           const VectorX<double>&,
                                           const VectorX<double>&,
                                           const VectorX<double>&,
                                           const VectorX<double>&);
template void normalize<double>(VectorX<double>&, VectorX<double>&,
                                VectorX<double>&);
template void ncalc<double>(const VectorX<double>&, const VectorX<double>&,
                            const VectorX<double>&, VectorX<double>&,
                            VectorX<double>&);
template TEParams<double> tecalc<double>(const VectorX<double>&,
                                         const VectorX<double>&,
                                         const VectorX<double>&,
                                         const VectorX<double>&, double);

}  // namespace xfoil
