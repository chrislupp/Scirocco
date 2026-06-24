/**
 * @brief Internal access to the per-closure _dRet/_dMa helpers and ThermoMaDeriv
 * for diagnostic probes (probe_per_closure_dma etc.). NOT a public API —
 * these change with the residual chain.
 */

#pragma once

#include "scirocco/config.h"
#include "scirocco/data/param.h"

namespace scirocco {
namespace detail {

struct ThermoMaDeriv {
    double KTb_dMa = 0;
    double KTl_dMa = 0;
    double H0_dMa = 0;
    double mu0_dMa = 0;
    double rho0_dMa = 0;
};

ThermoMaDeriv ComputeThermoMaDeriv(const Param& param, double Vinf,
                                   double rhoinf);

double GetUk_dMa(double ue, const Param& param, const ThermoMaDeriv& d);
double GetMach2_dMa(const Vector4d& U, const Param& param,
                    const ThermoMaDeriv& d);
double GetHk_dMa(const Vector4d& U, const Param& param, const ThermoMaDeriv& d);
double GetHss_dMa(const Vector4d& U, const Param& param,
                  const ThermoMaDeriv& d);
double GetRet_dMa(const Vector4d& U, const Param& param,
                  const ThermoMaDeriv& d);
double GetHs_dMa(const Vector4d& U, const Param& param, const ThermoMaDeriv& d);
double GetUs_dMa(const Vector4d& U, const Param& param, const ThermoMaDeriv& d);
double GetCf_dMa(const Vector4d& U, const Param& param, const ThermoMaDeriv& d);
double GetCDi_dMa(const Vector4d& U, const Param& param,
                  const ThermoMaDeriv& d);
double GetCfxt_dMa(const Vector4d& U, double x, const Param& param,
                   const ThermoMaDeriv& d);
double GetCDixt_dMa(const Vector4d& U, double x, const Param& param,
                    const ThermoMaDeriv& d);

}  // namespace detail
}  // namespace scirocco
