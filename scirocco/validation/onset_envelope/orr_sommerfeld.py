"""Orr-Sommerfeld spatial stability of Falkner-Skan profiles.

Lengths are scaled by the displacement thickness delta*, velocities by U_e,
and Re = U_e delta* / nu. Chebyshev collocation on an algebraically mapped
wall-normal grid; phi = phi' = 0 at the wall and the far field.

Checks (run this file): Jordinson's Blasius point Re = 998, alpha = 0.308,
c = 0.36412 + 0.00796i; Blasius critical Re = 519.4.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.linalg import eig

def falkner_skan(beta, eta_max=12.0):
    """f''' + f f'' + beta (1 - f'^2) = 0. Returns callables U(y), U''(y) in
    delta* units, plus H and theta_eta."""
    def rhs(e, z): return [z[1], z[2], -z[0]*z[2] - beta*(1 - z[1]**2)]
    def miss(s):
        sol = solve_ivp(rhs, [0, eta_max], [0, 0, s], rtol=1e-11, atol=1e-12)
        return sol.y[1, -1] - 1.0
    # bracket f''(0)
    grid = np.linspace(0.01 if beta < 0 else 0.2, 2.5, 120)
    vals = [miss(g) for g in grid]
    for a, b, va, vb in zip(grid[:-1], grid[1:], vals[:-1], vals[1:]):
        if np.isfinite(va) and np.isfinite(vb) and va*vb < 0:
            s = brentq(miss, a, b, xtol=1e-14); break
    else:
        raise RuntimeError('no bracket')
    sol = solve_ivp(rhs, [0, eta_max], [0, 0, s], rtol=1e-12, atol=1e-13, dense_output=True)
    e = np.linspace(0, eta_max, 6001); z = sol.sol(e)
    fp = z[1]
    dstar = np.trapezoid(1 - fp, e); theta = np.trapezoid(fp*(1 - fp), e)
    def U(y):
        ee = np.clip(np.asarray(y)*dstar, 0, eta_max); return sol.sol(ee)[1]
    def Upp(y):
        ee = np.clip(np.asarray(y)*dstar, 0, eta_max); zz = sol.sol(ee)
        fppp = -zz[0]*zz[2] - beta*(1 - zz[1]**2)
        return fppp*dstar**2
    return U, Upp, dstar/theta, theta, s

def cheb(N):
    x = np.cos(np.pi*np.arange(N+1)/N)
    c = np.r_[2, np.ones(N-1), 2]*(-1)**np.arange(N+1)
    X = np.tile(x, (N+1, 1)).T
    dX = X - X.T
    D = np.outer(c, 1/c)/(dX + np.eye(N+1))
    D -= np.diag(D.sum(axis=1))
    return D, x

class OS:
    def __init__(self, beta, N=90, ymax=40.0, yi=3.0):
        self.U, self.Upp, self.H, self.theta_eta, self.fpp0 = falkner_skan(beta)
        D, xi = cheb(N)
        # algebraic map y in [0, ymax] clustering half the points below yi
        a = yi*ymax/(ymax - 2*yi); b = 1 + 2*a/ymax
        y = a*(1 + xi)/(b - xi)
        dydxi = a*(b + 1)/(b - xi)**2
        D1 = np.diag(1/dydxi) @ D
        self.y = y; self.D1 = D1; self.D2 = D1 @ D1; self.D4 = self.D2 @ self.D2
        self.Uy = self.U(y); self.Uppy = self.Upp(y); self.N = N
        self.wall = int(np.argmin(y)); self.far = int(np.argmax(y))

    def omegas(self, alpha, Re):
        n = self.N + 1; I = np.eye(n)
        L = self.D2 - alpha**2*I
        A = alpha*np.diag(self.Uy) @ L - alpha*np.diag(self.Uppy) + (1j/Re)*(L @ L)
        B = L.astype(complex)
        # BCs phi = phi' = 0 at wall and far field: replace 4 rows
        # use rows adjacent to boundaries for the derivative conditions
        rw2 = self.wall - 1 if self.wall > 0 else self.wall + 1
        rf2 = self.far + 1 if self.far < self.N else self.far - 1
        A[self.wall] = I[self.wall]; B[self.wall] = 0
        A[self.far] = I[self.far]; B[self.far] = 0
        A[rw2] = self.D1[self.wall]; B[rw2] = 0
        A[rf2] = self.D1[self.far]; B[rf2] = 0
        w = eig(A, B, right=False)
        w = w[np.isfinite(w)]
        c = w/alpha
        ok = (c.real > 0.05) & (c.real < 0.95)
        return w[ok]

    def ts_omega(self, alpha, Re, guess=None):
        w = self.omegas(alpha, Re)
        if guess is None:
            return w[np.argmax(w.imag)]
        return w[np.argmin(np.abs(w - guess))]

    def spatial_alpha(self, omega, Re, alpha0, w_guess=None, tol=1e-10):
        """Find complex alpha with TS omega(alpha) == omega (real)."""
        a0 = alpha0; f0 = self.ts_omega(a0, Re, w_guess) - omega
        a1 = alpha0*(1 + 1e-3) + 1e-4j; f1 = self.ts_omega(a1, Re, f0 + omega) - omega
        for _ in range(40):
            if abs(f1 - f0) == 0: break
            a2 = a1 - f1*(a1 - a0)/(f1 - f0)
            a0, f0 = a1, f1
            a1 = a2; f1 = self.ts_omega(a1, Re, f0 + omega) - omega
            if abs(f1) < tol: return a1
        if abs(f1) < 1e-7: return a1
        raise RuntimeError('spatial solve failed omega=%g Re=%g' % (omega, Re))


if __name__ == '__main__':
    from scipy.optimize import minimize_scalar
    s = OS(0.0)
    print('Blasius H = %.4f (ref 2.5911)' % s.H)
    print('c(Re=998, alpha=0.308) =', s.ts_omega(0.308, 998.0)/0.308, '(ref 0.36412+0.00796j)')
    def wi_max(Re):
        r = minimize_scalar(lambda a: -s.ts_omega(a, Re).imag, bounds=(0.15, 0.4), method='bounded', options={'xatol': 1e-7})
        return -r.fun
    print('Re_crit = %.1f (ref 519.4)' % brentq(wi_max, 400, 600, xtol=1e-4))
