---
sidebar_position: 7
---

# `KulfanDiscipline`

A Philote-Cpp explicit discipline that maps Kulfan / CST Bernstein
coefficients to panel-node coordinates with analytic partials. Used
inside MDO frameworks (OpenMDAO via `philote-mdo`) to expose CST
parameterization as an upstream component to the Scirocco analysis
discipline. Header:
[`scirocco/kulfan_discipline.h`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/include/scirocco/kulfan_discipline.h).

This is also the discipline implemented by the `scirocco_kulfan_server`
binary — see [Optimization → Server](../optimization/server) for the
runtime story.

## Class shape

Each surface is parameterized as

$$
y_\text{upper}(\psi) = +\,C(\psi)\,S_\text{upper}(\psi) + \psi\, y_{\text{TE},\text{upper}}
$$

$$
y_\text{lower}(\psi) = -\,C(\psi)\,S_\text{lower}(\psi) + \psi\, y_{\text{TE},\text{lower}}
$$

with $\psi = x/c \in [0, 1]$, class function

$$
C(\psi) = \psi^{1/2} (1-\psi)^{1}
$$

(round LE, sharp TE), and shape function

$$
S(\psi) = \sum_{i=0}^{n} A_i \binom{n}{i} \psi^i (1-\psi)^{n-i}.
$$

The Bernstein coefficients $A_i$ are the design variables. With the sign
convention above, $A_i \ge 0$ on both surfaces gives a non-self-
intersecting airfoil **everywhere** on $[0, 1]$, not just at sample
points.

## Options

Set at discipline construction.

| Option | Type | Default | Meaning |
| --- | --- | --- | --- |
| `n_cst_upper` | int | 6 | Number of upper-surface Bernstein coefficients |
| `n_cst_lower` | int | 6 | Number of lower-surface Bernstein coefficients |
| `npanel` | int (even) | 200 | Total panels; `npoint = npanel + 1` |
| `te_upper` | double | 0.0 | Fixed upper-surface TE y-offset |
| `te_lower` | double | 0.0 | Fixed lower-surface TE y-offset |

## Inputs

| Name | Shape |
| --- | --- |
| `cst_upper` | `n_cst_upper` |
| `cst_lower` | `n_cst_lower` |

## Outputs

| Name | Shape | Meaning |
| --- | --- | --- |
| `x` | `2 · npoint` | Panel-node coordinates, x-block then z-block, TE-lower → LE → TE-upper. **Directly compatible with `SiroccoDiscipline.x` input** |
| `smoothness_upper` | `n_cst_upper − 2` | Second differences $D_i = A_{i+2} - 2A_{i+1} + A_i$ for the upper-surface coefficients |
| `smoothness_lower` | `n_cst_lower − 2` | Same, lower surface |

The smoothness outputs are the right knob for stabilizing a CST viscous
optimizer — bounding $|D_i|$ caps the second derivative of $S(\psi)$, and
hence the airfoil-surface curvature, with a linear constraint. See
[CST design](../optimization/cst-design) for the recipe.

## C++ class

```cpp
class KulfanDiscipline : public philote::ExplicitDiscipline {
public:
    KulfanDiscipline() = default;
    ~KulfanDiscipline() noexcept = default;

    void Initialize()                                                  override;
    void SetOptions(const google::protobuf::Struct& options)           override;
    void Configure()                                                   override;
    void Setup()                                                       override;
    void SetupPartials()                                               override;
    void Compute(const philote::Variables& inputs,
                 philote::Variables&        outputs)                    override;
    void ComputePartials(const philote::Variables& inputs,
                         philote::Partials&         partials)           override;

    // Test-friendly accessors
    int               npoint()      const;
    int               npanel()      const;
    int               n_cst_upper() const;
    int               n_cst_lower() const;
    const VectorXd&   psi_eval()    const;
};
```

You typically don't instantiate `KulfanDiscipline` directly from C++ —
you launch the `scirocco_kulfan_server` binary and talk to it via gRPC
from Python (`RemoteExplicitComponent`).

## Evaluation stations

The discipline samples the surface at cosine-clustered $\psi$ stations
on the half-chord. The number of stations is `npanel/2 + 1`, with
mirror points giving the full panel-node grid. The `psi_eval()` accessor
exposes the half-chord $\psi$ vector for testing.

## Building, launching, using

See [Optimization → Server](../optimization/server) for the build
recipe and gRPC launch commands. Once the server is up, use it from
OpenMDAO as a `RemoteExplicitComponent`:

```python
from philote_mdo.openmdao import RemoteExplicitComponent
import grpc

channel = grpc.insecure_channel("localhost:50052")
geometry = RemoteExplicitComponent(
    channel=channel,
    n_cst_upper=8, n_cst_lower=8, npanel=100,
    te_upper=+0.003, te_lower=-0.003,
)
```

## See also

- Tutorial: not yet — the canonical CST tutorial is the inviscid + viscous
  optimization walkthrough in [CST design](../optimization/cst-design)
