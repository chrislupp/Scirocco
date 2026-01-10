// Copyright (c) 2024 ProjectTrysil Authors
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Airfoil geometry utilities.

#include "xfoil/geometry.h"

#include <algorithm>
#include <cmath>

#include "xfoil/spline.h"

namespace xfoil {

NacaResult naca4(int ides, int nside) {
    // Parse digits
    int n4 = ides / 1000;
    int n3 = (ides - n4 * 1000) / 100;
    int n2 = (ides - n4 * 1000 - n3 * 100) / 10;
    int n1 = ides - n4 * 1000 - n3 * 100 - n2 * 10;

    double m = static_cast<double>(n4) / 100.0;
    double p = static_cast<double>(n3) / 10.0;
    double t = static_cast<double>(n2 * 10 + n1) / 100.0;

    NacaResult result = naca4_params(m, p, t, nside);
    result.name = "NACA " + std::to_string(n4) + std::to_string(n3) +
                  std::to_string(n2) + std::to_string(n1);
    return result;
}

NacaResult naca4_params(double m, double p, double t, int nside) {
    NacaResult result;

    // TE point bunching parameter
    const double an = 1.5;
    const double anp = an + 1.0;

    VectorXd xx(nside), yt(nside), yc(nside);

    for (int i = 0; i < nside; ++i) {
        double frac = static_cast<double>(i) / static_cast<double>(nside - 1);
        if (i == nside - 1) {
            xx(i) = 1.0;
        } else {
            xx(i) = 1.0 - anp * frac * std::pow(1.0 - frac, an) -
                    std::pow(1.0 - frac, anp);
        }

        // Thickness distribution
        yt(i) = (0.29690 * std::sqrt(xx(i)) - 0.12600 * xx(i) -
                 0.35160 * xx(i) * xx(i) +
                 0.28430 * xx(i) * xx(i) * xx(i) -
                 0.10150 * xx(i) * xx(i) * xx(i) * xx(i)) *
                t / 0.20;

        // Camber distribution
        if (p > 0.0) {
            if (xx(i) < p) {
                yc(i) = m / (p * p) * (2.0 * p * xx(i) - xx(i) * xx(i));
            } else {
                yc(i) = m / ((1.0 - p) * (1.0 - p)) *
                        ((1.0 - 2.0 * p) + 2.0 * p * xx(i) - xx(i) * xx(i));
            }
        } else {
            yc(i) = 0.0;
        }
    }

    // Build airfoil coordinates
    int nb = 2 * nside - 1;
    result.x.resize(nb);
    result.y.resize(nb);

    int ib = 0;
    for (int i = nside - 1; i >= 0; --i) {
        result.x(ib) = xx(i);
        result.y(ib) = yc(i) + yt(i);
        ++ib;
    }
    for (int i = 1; i < nside; ++i) {
        result.x(ib) = xx(i);
        result.y(ib) = yc(i) - yt(i);
        ++ib;
    }

    result.n = nb;
    result.name = "NACA airfoil";
    return result;
}

NacaResult naca5(int ides, int nside) {
    NacaResult result;

    // Parse digits
    int n5 = ides / 10000;
    int n4 = (ides - n5 * 10000) / 1000;
    int n3 = (ides - n5 * 10000 - n4 * 1000) / 100;
    int n2 = (ides - n5 * 10000 - n4 * 1000 - n3 * 100) / 10;
    int n1 = ides - n5 * 10000 - n4 * 1000 - n3 * 100 - n2 * 10;

    int n543 = 100 * n5 + 10 * n4 + n3;

    // Look up camber line parameters
    double m_cl, c_cl;
    if (n543 == 210) {
        m_cl = 0.0580;
        c_cl = 361.4;
    } else if (n543 == 220) {
        m_cl = 0.1260;
        c_cl = 51.64;
    } else if (n543 == 230) {
        m_cl = 0.2025;
        c_cl = 15.957;
    } else if (n543 == 240) {
        m_cl = 0.2900;
        c_cl = 6.643;
    } else if (n543 == 250) {
        m_cl = 0.3910;
        c_cl = 3.230;
    } else {
        result.n = 0;
        return result;
    }

    double t = static_cast<double>(n2 * 10 + n1) / 100.0;

    const double an = 1.5;
    const double anp = an + 1.0;

    VectorXd xx(nside), yt(nside), yc(nside);

    for (int i = 0; i < nside; ++i) {
        double frac = static_cast<double>(i) / static_cast<double>(nside - 1);
        if (i == nside - 1) {
            xx(i) = 1.0;
        } else {
            xx(i) = 1.0 - anp * frac * std::pow(1.0 - frac, an) -
                    std::pow(1.0 - frac, anp);
        }

        yt(i) = (0.29690 * std::sqrt(xx(i)) - 0.12600 * xx(i) -
                 0.35160 * xx(i) * xx(i) +
                 0.28430 * xx(i) * xx(i) * xx(i) -
                 0.10150 * xx(i) * xx(i) * xx(i) * xx(i)) *
                t / 0.20;

        if (xx(i) < m_cl) {
            yc(i) = (c_cl / 6.0) *
                    (xx(i) * xx(i) * xx(i) -
                     3.0 * m_cl * xx(i) * xx(i) +
                     m_cl * m_cl * (3.0 - m_cl) * xx(i));
        } else {
            yc(i) = (c_cl / 6.0) * m_cl * m_cl * m_cl * (1.0 - xx(i));
        }
    }

    int nb = 2 * nside - 1;
    result.x.resize(nb);
    result.y.resize(nb);

    int ib = 0;
    for (int i = nside - 1; i >= 0; --i) {
        result.x(ib) = xx(i);
        result.y(ib) = yc(i) + yt(i);
        ++ib;
    }
    for (int i = 1; i < nside; ++i) {
        result.x(ib) = xx(i);
        result.y(ib) = yc(i) - yt(i);
        ++ib;
    }

    result.n = nb;
    result.name = "NACA " + std::to_string(n5) + std::to_string(n4) +
                  std::to_string(n3) + std::to_string(n2) + std::to_string(n1);
    return result;
}

double lefind(const VectorXd& x, const VectorXd& xp,
              const VectorXd& y, const VectorXd& yp,
              const VectorXd& s) {
    int n = x.size();

    double dseps = (s(n - 1) - s(0)) * 1.0e-5;

    double xte = 0.5 * (x(0) + x(n - 1));
    double yte = 0.5 * (y(0) + y(n - 1));

    int ile = 2;
    for (int i = 2; i < n - 2; ++i) {
        double dxte = x(i) - xte;
        double dyte = y(i) - yte;
        double dx = x(i + 1) - x(i);
        double dy = y(i + 1) - y(i);
        double dotp = dxte * dx + dyte * dy;
        if (dotp < 0.0) {
            ile = i;
            break;
        }
    }

    double sle = s(ile);

    if (s(ile) == s(ile - 1)) {
        return sle;
    }

    // Newton iteration
    for (int iter = 0; iter < 50; ++iter) {
        double xle = seval(sle, x, xp, s);
        double yle = seval(sle, y, yp, s);
        double dxds = deval(sle, x, xp, s);
        double dyds = deval(sle, y, yp, s);
        double dxdd = d2val(sle, x, xp, s);
        double dydd = d2val(sle, y, yp, s);

        double xchord = xle - xte;
        double ychord = yle - yte;

        double res = xchord * dxds + ychord * dyds;
        double ress = dxds * dxds + dyds * dyds + xchord * dxdd + ychord * dydd;

        double dsle = -res / ress;

        double chord_approx = std::abs(xchord + ychord);
        dsle = std::max(dsle, -0.02 * chord_approx);
        dsle = std::min(dsle, 0.02 * chord_approx);
        sle = sle + dsle;

        if (std::abs(dsle) < dseps) {
            return sle;
        }
    }

    return s(ile);
}

GeomParams geopar(const VectorXd& x, const VectorXd& xp,
                  const VectorXd& y, const VectorXd& yp,
                  const VectorXd& s) {
    GeomParams result;
    int n = x.size();

    result.sle = lefind(x, xp, y, yp, s);
    result.xle = seval(result.sle, x, xp, s);
    result.yle = seval(result.sle, y, yp, s);
    result.xte = 0.5 * (x(0) + x(n - 1));
    result.yte = 0.5 * (y(0) + y(n - 1));

    double dx = result.xte - result.xle;
    double dy = result.yte - result.yle;
    result.chord = std::sqrt(dx * dx + dy * dy);

    double curvle = curv(result.sle, x, xp, y, yp, s);
    result.radle = 0.0;
    if (std::abs(curvle) > 0.001 * (s(n - 1) - s(0))) {
        result.radle = 1.0 / curvle;
    }

    double ang1 = std::atan2(-yp(0), -xp(0));
    double ang2 = std::atan2(yp(n - 1), xp(n - 1));
    result.angte = ang2 - ang1;

    return result;
}

void normalize(VectorXd& x, VectorXd& y, VectorXd& s) {
    int n = x.size();

    scalc(x, y, s);
    VectorXd xp, yp;
    segspl(x, s, xp);
    segspl(y, s, yp);

    double sle = lefind(x, xp, y, yp, s);

    double xmin = seval(sle, x, xp, s);
    double ymin = seval(sle, y, yp, s);
    double xmax = 0.5 * (x(0) + x(n - 1));

    double fudge = 1.0 / (xmax - xmin);

    for (int i = 0; i < n; ++i) {
        x(i) = (x(i) - xmin) * fudge;
        y(i) = (y(i) - ymin) * fudge;
        s(i) = s(i) * fudge;
    }
}

void ncalc(const VectorXd& x, const VectorXd& y, const VectorXd& s,
           VectorXd& nx, VectorXd& ny) {
    int n = x.size();
    if (n <= 1) return;

    nx.resize(n);
    ny.resize(n);

    segspl(x, s, nx);
    segspl(y, s, ny);

    for (int i = 0; i < n; ++i) {
        double sx = ny(i);
        double sy = -nx(i);
        double smod = std::sqrt(sx * sx + sy * sy);
        nx(i) = sx / smod;
        ny(i) = sy / smod;
    }

    for (int i = 0; i < n - 1; ++i) {
        if (s(i) == s(i + 1)) {
            double sx = 0.5 * (nx(i) + nx(i + 1));
            double sy = 0.5 * (ny(i) + ny(i + 1));
            double smod = std::sqrt(sx * sx + sy * sy);
            nx(i) = sx / smod;
            ny(i) = sy / smod;
            nx(i + 1) = nx(i);
            ny(i + 1) = ny(i);
        }
    }
}

TEParams tecalc(const VectorXd& x, const VectorXd& xp,
                const VectorXd& y, const VectorXd& yp, double chord) {
    TEParams te;
    int n = x.size();

    te.dxte = x(0) - x(n - 1);
    te.dyte = y(0) - y(n - 1);

    double dxs = 0.5 * (-xp(0) + xp(n - 1));
    double dys = 0.5 * (-yp(0) + yp(n - 1));

    te.ante = dxs * te.dyte - dys * te.dxte;
    te.aste = dxs * te.dxte + dys * te.dyte;

    te.dste = std::sqrt(te.dxte * te.dxte + te.dyte * te.dyte);
    te.sharp = te.dste < 0.0001 * chord;

    if (te.sharp) {
        te.scs = 1.0;
        te.sds = 0.0;
    } else {
        te.scs = te.ante / te.dste;
        te.sds = te.aste / te.dste;
    }

    return te;
}

PanelResult pangen(const VectorXd& xb, const VectorXd& yb,
                   const PanelParams& params) {
    PanelResult result;
    int nb = xb.size();

    if (nb < 2) {
        result.n = 0;
        result.converged = false;
        return result;
    }

    const int ipfac = 3;
    int n = params.npan;

    VectorXd sb, xbp, ybp;
    scalc(xb, yb, sb);
    segspl(xb, sb, xbp);
    segspl(yb, sb, ybp);

    double sbref = 0.5 * (sb(nb - 1) - sb(0));

    VectorXd w5(nb);
    for (int i = 0; i < nb; ++i) {
        w5(i) = std::abs(curv(sb(i), xb, xbp, yb, ybp, sb)) * sbref;
    }

    double sble = lefind(xb, xbp, yb, ybp, sb);
    double cvle = std::abs(curv(sble, xb, xbp, yb, ybp, sb)) * sbref;

    int ible = 0;
    for (int i = 0; i < nb - 1; ++i) {
        if (sble == sb(i) && sble == sb(i + 1)) {
            ible = i + 1;
            break;
        }
    }

    double xble = seval(sble, xb, xbp, sb);
    double yble = seval(sble, yb, ybp, sb);
    double xbte = 0.5 * (xb(0) + xb(nb - 1));
    double ybte = 0.5 * (yb(0) + yb(nb - 1));
    double chbsq = (xbte - xble) * (xbte - xble) + (ybte - yble) * (ybte - yble);

    int nk = 3;
    double cvsum = 0.0;
    for (int k = -nk; k <= nk; ++k) {
        double frac = static_cast<double>(k) / static_cast<double>(nk);
        double sbk = sble + frac * sbref / std::max(cvle, 20.0);
        double cvk = std::abs(curv(sbk, xb, xbp, yb, ybp, sb)) * sbref;
        cvsum += cvk;
    }
    double cvavg = cvsum / static_cast<double>(2 * nk + 1);
    if (ible != 0) cvavg = 10.0;

    double cc = 6.0 * params.cvpar;
    double cvte = cvavg * params.cterat;
    w5(0) = cvte;
    w5(nb - 1) = cvte;

    double smool = std::max(1.0 / std::max(cvavg, 20.0),
                            0.25 / static_cast<double>(n / 2));
    double smoosq = (smool * sbref) * (smool * sbref);

    VectorXd w1(nb), w2(nb), w3(nb);
    w2(0) = 1.0;
    w3(0) = 0.0;

    for (int i = 1; i < nb - 1; ++i) {
        double dsm = sb(i) - sb(i - 1);
        double dsp = sb(i + 1) - sb(i);
        double dso = 0.5 * (sb(i + 1) - sb(i - 1));

        if (dsm == 0.0 || dsp == 0.0) {
            w1(i) = 0.0;
            w2(i) = 1.0;
            w3(i) = 0.0;
        } else {
            w1(i) = smoosq * (-1.0 / dsm) / dso;
            w2(i) = smoosq * (1.0 / dsp + 1.0 / dsm) / dso + 1.0;
            w3(i) = smoosq * (-1.0 / dsp) / dso;
        }
    }
    w1(nb - 1) = 0.0;
    w2(nb - 1) = 1.0;

    for (int i = 1; i < nb - 1; ++i) {
        if (sb(i) == sble || i == ible - 1 || i == ible) {
            w1(i) = 0.0;
            w2(i) = 1.0;
            w3(i) = 0.0;
            w5(i) = cvle;
        } else if (sb(i - 1) < sble && sb(i) > sble) {
            double dsm = sb(i - 1) - sb(i - 2);
            double dsp = sble - sb(i - 1);
            double dso = 0.5 * (sble - sb(i - 2));

            w1(i - 1) = smoosq * (-1.0 / dsm) / dso;
            w2(i - 1) = smoosq * (1.0 / dsp + 1.0 / dsm) / dso + 1.0;
            w3(i - 1) = 0.0;
            w5(i - 1) = w5(i - 1) + smoosq * cvle / (dsp * dso);

            dsm = sb(i) - sble;
            dsp = sb(i + 1) - sb(i);
            dso = 0.5 * (sb(i + 1) - sble);
            w1(i) = 0.0;
            w2(i) = smoosq * (1.0 / dsp + 1.0 / dsm) / dso + 1.0;
            w3(i) = smoosq * (-1.0 / dsp) / dso;
            w5(i) = w5(i) + smoosq * cvle / (dsm * dso);
            break;
        }
    }

    for (int i = 1; i < nb - 1; ++i) {
        double xoc = ((xb(i) - xble) * (xbte - xble) +
                      (yb(i) - yble) * (ybte - yble)) / chbsq;
        if (sb(i) < sble) {
            if (xoc > params.xsref1 && xoc < params.xsref2) {
                w1(i) = 0.0;
                w2(i) = 1.0;
                w3(i) = 0.0;
                w5(i) = cvle * params.ctrrat;
            }
        } else {
            if (xoc > params.xpref1 && xoc < params.xpref2) {
                w1(i) = 0.0;
                w2(i) = 1.0;
                w3(i) = 0.0;
                w5(i) = cvle * params.ctrrat;
            }
        }
    }

    if (ible == 0) {
        trisol(w2, w1, w3, w5);
    } else {
        VectorXd w2a = w2.head(ible);
        VectorXd w1a = w1.head(ible);
        VectorXd w3a = w3.head(ible);
        VectorXd w5a = w5.head(ible);
        trisol(w2a, w1a, w3a, w5a);
        w5.head(ible) = w5a;

        VectorXd w2b = w2.tail(nb - ible);
        VectorXd w1b = w1.tail(nb - ible);
        VectorXd w3b = w3.tail(nb - ible);
        VectorXd w5b = w5.tail(nb - ible);
        trisol(w2b, w1b, w3b, w5b);
        w5.tail(nb - ible) = w5b;
    }

    double cvmax = w5.cwiseAbs().maxCoeff();
    w5 /= cvmax;

    VectorXd w6;
    segspl(w5, sb, w6);

    int nn = ipfac * (n - 1) + 1;
    VectorXd snew(nn);

    double rdste = 0.667;
    double rtf = (rdste - 1.0) * static_cast<double>(ipfac) + 1.0;

    if (ible == 0) {
        double dsavg = (sb(nb - 1) - sb(0)) /
                       (static_cast<double>(nn - 3) + 2.0 * rtf);
        snew(0) = sb(0);
        for (int i = 1; i < nn - 1; ++i) {
            snew(i) = sb(0) + dsavg * (static_cast<double>(i - 1) + rtf);
        }
        snew(nn - 1) = sb(nb - 1);
    } else {
        int nfrac1 = (n * ible) / nb;
        int nn1 = ipfac * (nfrac1 - 1) + 1;
        double dsavg1 = (sble - sb(0)) / (static_cast<double>(nn1 - 2) + rtf);
        snew(0) = sb(0);
        for (int i = 1; i < nn1; ++i) {
            snew(i) = sb(0) + dsavg1 * (static_cast<double>(i - 1) + rtf);
        }

        int nn2 = nn - nn1 + 1;
        double dsavg2 = (sb(nb - 1) - sble) / (static_cast<double>(nn2 - 2) + rtf);
        for (int i = 1; i < nn2 - 1; ++i) {
            snew(i - 1 + nn1) = sble + dsavg2 * (static_cast<double>(i - 1) + rtf);
        }
        snew(nn - 1) = sb(nb - 1);
    }

    VectorXd w4(nn);
    w1.resize(nn);
    w2.resize(nn);
    w3.resize(nn);
    result.converged = false;

    for (int iter = 0; iter < 20; ++iter) {
        double cv1 = seval(snew(0), w5, w6, sb);
        double cv2 = seval(snew(1), w5, w6, sb);
        double cvs1 = deval(snew(0), w5, w6, sb);
        double cvs2 = deval(snew(1), w5, w6, sb);

        double cavm = std::sqrt(cv1 * cv1 + cv2 * cv2);
        double cavm_s1 = (cavm == 0.0) ? 0.0 : cvs1 * cv1 / cavm;
        double cavm_s2 = (cavm == 0.0) ? 0.0 : cvs2 * cv2 / cavm;

        for (int i = 1; i < nn - 1; ++i) {
            double dsm = snew(i) - snew(i - 1);
            double dsp = snew(i) - snew(i + 1);
            double cv3 = seval(snew(i + 1), w5, w6, sb);
            double cvs3 = deval(snew(i + 1), w5, w6, sb);

            double cavp = std::sqrt(cv3 * cv3 + cv2 * cv2);
            double cavp_s2 = (cavp == 0.0) ? 0.0 : cvs2 * cv2 / cavp;
            double cavp_s3 = (cavp == 0.0) ? 0.0 : cvs3 * cv3 / cavp;

            double fm = cc * cavm + 1.0;
            double fp = cc * cavp + 1.0;

            double rez = dsp * fp + dsm * fm;

            w1(i) = -fm + cc * dsm * cavm_s1;
            w2(i) = fp + fm + cc * (dsp * cavp_s2 + dsm * cavm_s2);
            w3(i) = -fp + cc * dsp * cavp_s3;
            w4(i) = -rez;

            cv1 = cv2;
            cv2 = cv3;
            cvs1 = cvs2;
            cvs2 = cvs3;
            cavm = cavp;
            cavm_s1 = cavp_s2;
            cavm_s2 = cavp_s3;
        }

        w2(0) = 1.0;
        w3(0) = 0.0;
        w4(0) = 0.0;
        w1(nn - 1) = 0.0;
        w2(nn - 1) = 1.0;
        w4(nn - 1) = 0.0;

        if (rtf != 1.0) {
            int i = 1;
            w4(i) = -((snew(i) - snew(i - 1)) + rtf * (snew(i) - snew(i + 1)));
            w1(i) = -1.0;
            w2(i) = 1.0 + rtf;
            w3(i) = -rtf;

            i = nn - 2;
            w4(i) = -((snew(i) - snew(i + 1)) + rtf * (snew(i) - snew(i - 1)));
            w3(i) = -1.0;
            w2(i) = 1.0 + rtf;
            w1(i) = -rtf;
        }

        trisol(w2, w1, w3, w4);

        double rlx = 1.0;
        double dmax = 0.0;
        for (int i = 0; i < nn - 1; ++i) {
            double ds = snew(i + 1) - snew(i);
            double dds = w4(i + 1) - w4(i);
            double dsrat = 1.0 + rlx * dds / ds;
            if (dsrat > 4.0) rlx = (4.0 - 1.0) * ds / dds;
            if (dsrat < 0.2) rlx = (0.2 - 1.0) * ds / dds;
            dmax = std::max(std::abs(w4(i)), dmax);
        }

        for (int i = 1; i < nn - 1; ++i) {
            snew(i) = snew(i) + rlx * w4(i);
        }

        if (std::abs(dmax) < 1.0e-3) {
            result.converged = true;
            break;
        }
    }

    result.x.resize(n);
    result.y.resize(n);
    result.s.resize(n);

    for (int i = 0; i < n; ++i) {
        int ind = ipfac * i;
        result.s(i) = snew(ind);
        result.x(i) = seval(snew(ind), xb, xbp, sb);
        result.y(i) = seval(snew(ind), yb, ybp, sb);
    }
    result.n = n;

    scalc(result.x, result.y, result.s);

    VectorXd xp, yp;
    segspl(result.x, result.s, xp);
    segspl(result.y, result.s, yp);
    result.sle = lefind(result.x, xp, result.y, yp, result.s);

    return result;
}

}  // namespace xfoil
