# Verification and validation data

Scripts behind the documentation pages
[XFOIL and mfoil](../../../docs/docs/validation/xfoil-mfoil.md) and
[Wind-tunnel data](../../../docs/docs/validation/experiments.md).

## Reproducing

Build the library in Release, then the driver against it:

```bash
cmake -S scirocco -B scirocco/build-release -DCMAKE_BUILD_TYPE=Release
cmake --build scirocco/build-release -j --target scirocco
c++ -O3 -std=c++17 -I scirocco/include -I scirocco/build-release/include \
    -isystem /opt/homebrew/include/eigen3 scirocco/validation/vv/vv_polar.cpp \
    scirocco/build-release/libscirocco.a -framework Accelerate -o vv_polar
```

On Linux, replace `-framework Accelerate` with the BLAS and LAPACK libraries
the library was configured with (for example `-lopenblas`), or configure
with `-DSCIROCCO_USE_BLAS=OFF`.

Then, with XFOIL 6.99 on the `PATH` (or set `XFOIL`) and NumPy and SciPy
installed for mfoil:

```bash
python3 scirocco/validation/vv/vv.py codes ./vv_polar docs/src/data/vv_codes.json
python3 scirocco/validation/vv/vv.py experiments ./vv_polar docs/src/data/vv_experiments.json
```

`vv.py` holds the case lists. `vv_polar.cpp` documents the driver's
arguments.

## Measurements

`experiments/` holds the measured data, unchanged, as distributed by the
[NASA Turbulence Modeling Resource](https://tmbwg.github.io/turbmodels/):

| File | Source |
| --- | --- |
| `0012.abbottdata.cl.dat`, `0012.abbottdata.cd.dat` | Abbott & von Doenhoff, *Theory of Wing Sections*, digitized |
| `CLCD_Ladson_expdata.dat` | Ladson, NASA TM-4074, 1988 |
| `CL_Gregory_expdata.dat`, `CP_Gregory_expdata.dat` | Gregory & O'Reilly, ARC R&M 3726, 1970, digitized |
| `naca4412.cp.expt.dat` | Coles & Wadcock, AIAA J. 17(4), 1979 |
