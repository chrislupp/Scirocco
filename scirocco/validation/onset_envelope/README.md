# Amplification onset from Orr–Sommerfeld envelopes

Data and fit behind `OnsetModel::Envelope` (C++) and `param.onset = 'envelope'`
(mfoil). The model is an opt-in alternative to the XFOIL/mfoil cubic onset
ramp.

## Why

The e^N model in XFOIL and mfoil grows the amplification factor as
`dn/dRe_theta = rf * da(Hk)`. `da` is the asymptotic envelope slope and `rf` is
an onset ramp. The default ramp is a cubic over `log10(Re_theta) = lrc(Hk) ± 0.1`.
It is a numerical smoothing device and is not fitted to stability data. It
switches growth on so sharply that shape gradients swing in sign from node to
node at the onset (issue #22), and refining the paneling does not remove the
swing.

## Method

`orr_sommerfeld.py` solves the spatial Orr–Sommerfeld problem for
Falkner–Skan profiles in displacement-thickness units. It uses Chebyshev
collocation on a mapped grid and finds complex α by a secant iteration on the
temporal eigenvalue. Run it to check the solver against Jordinson's Blasius
point (Re = 998, α = 0.308, c = 0.36412 + 0.00796i) and the Blasius critical
Reynolds number of 519.4.

`envelope.py BETA KMAX` computes the envelope N-factor curve against Re_theta
for one profile and writes `data/env_<beta>.json`. At each Re_δ* it anchors
on the temporally most unstable TS mode and sweeps ω out to both neutral
points. Each frequency is then integrated at constant
`F = ω ν / U_e²`, using `dN/dRe_theta = -α_i / (H θ_η²)` from Falkner–Skan
similarity. The stored runs cover β = −0.19 to 0.3, which is H = 2.36 to 3.48.
At β = 0.5 the curve never reaches the asymptotic range (N only reaches 4), so
that run is omitted.

`fit_onset.py` normalizes each envelope slope by its fitted asymptote and fits
`rf = P(3, k (z − zc))` against `z = log10(Re_theta / R0)`. P is the
regularized lower incomplete gamma function, so `rf` is exactly zero before
onset, C² at onset, and tends to 1. The rate `k = ln10 / (1 − 10^(zc/3))`
comes from `E[10^(zc+T)] = 1` with T ~ Gamma(3, k). That constraint keeps the
asymptote `n = da (Re_theta − 10^lrc)` identical to the cubic ramp's, so only
the onset changes. Pass `--docs PATH` to also write the documentation figure
data.

## Results

```
 beta     H    slope/da  log10R0-lrc  zcrit   zc_fit  rms_env  rms_cubic
+0.300  2.362   0.675      -0.028      -0.116  -0.159   0.056    0.126
+0.200  2.411   0.778      -0.051      -0.144  -0.202   0.047    0.153
+0.100  2.481   0.886      -0.032      -0.178  -0.260   0.034    0.174
+0.000  2.591   0.989      +0.040      -0.189  -0.282   0.020    0.184
-0.050  2.676   1.044      +0.062      -0.189  -0.276   0.018    0.180
-0.100  2.801   1.091      +0.052      -0.208  -0.300   0.018    0.192
-0.140  2.963   1.103      +0.025      -0.240  -0.343   0.019    0.213
-0.170  3.178   1.084      -0.001      -0.271  -0.380   0.022    0.229
-0.190  3.481   1.045      -0.017      -0.293  -0.404   0.025    0.237

zc = -0.6877 + 0.6797 / (Hk - 1)   (max residual 0.031)
```

- The computed envelopes support the correlation the model keeps. The
  asymptotic intercept agrees with `lrc` within 0.06 decades. The slope agrees
  with `da` within 10% for H ≥ 2.48. For stronger favorable gradients the
  slope is 10–30% below `da`. The onset model leaves `da` unchanged.
- Growth starts at the true critical Re_theta, 0.12–0.29 decades below the
  intercept. The cubic ramp starts only 0.1 decades below it.
- The envelope slope is asymmetric about the intercept: a quick rise to about
  0.65 at the intercept, then a long tail that reaches 1 near +0.4. The
  gamma-CDF ramp follows that shape with rms 0.02–0.06, compared with
  0.13–0.24 for the cubic.
- The fit applies the shape about `lrc`, not about each profile's own R0. A
  smooth softplus caps `zc` below −0.05 so that the constraint on `k` always
  has a solution.

## Effect on gradients

NACA 2412, α = 4°, Re = 10⁶. The values are dCd/dA for a Gaussian normal bump
of width 0.01c on the upper surface, swept through the onset:

| | peak positive | peak negative | secondary lobe |
| --- | --- | --- | --- |
| cubic (default) | +0.77 | −0.64 | +0.10 |
| envelope | +0.29 | −0.15 | none |

Both converge under refinement (160 to 480 panels). Cd changes by −0.06% and
upper-surface transition moves 0.0007c aft.

## Assumptions

- Locally parallel, two-dimensional TS waves; the most unstable mode only.
- Falkner–Skan similarity is used to map x to Re_theta.
- The envelope is taken over constant `F = ω ν / U_e²`, which is the usual
  e^N convention.
