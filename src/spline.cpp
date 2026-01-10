// Copyright (c) 2024 ProjectTrysil Authors
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Spline interpolation utilities - explicit template instantiations.

#include "xfoil/spline.h"

namespace xfoil {

// Explicit template instantiations for double
template void trisol<double>(VectorX<double>&, const VectorX<double>&,
                             VectorX<double>&, VectorX<double>&);
template void spline<double>(const VectorX<double>&, const VectorX<double>&,
                             VectorX<double>&);
template void splind<double>(const VectorX<double>&, const VectorX<double>&,
                             VectorX<double>&, double, double);
template void splina<double>(const VectorX<double>&, const VectorX<double>&,
                             VectorX<double>&);
template int find_interval<double>(double, const VectorX<double>&);
template double seval<double>(double, const VectorX<double>&,
                              const VectorX<double>&, const VectorX<double>&);
template double deval<double>(double, const VectorX<double>&,
                              const VectorX<double>&, const VectorX<double>&);
template double d2val<double>(double, const VectorX<double>&,
                              const VectorX<double>&, const VectorX<double>&);
template double curv<double>(double, const VectorX<double>&,
                             const VectorX<double>&, const VectorX<double>&,
                             const VectorX<double>&, const VectorX<double>&);
template bool sinvrt<double>(double&, double, const VectorX<double>&,
                             const VectorX<double>&, const VectorX<double>&);
template void scalc<double>(const VectorX<double>&, const VectorX<double>&,
                            VectorX<double>&);
template void segspl<double>(const VectorX<double>&, const VectorX<double>&,
                             VectorX<double>&);

}  // namespace xfoil
