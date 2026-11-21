---
sidebar_position: 7
---

# `KulfanDiscipline`

The Philote discipline that turns CST coefficients into panel nodes, with
exact partials. `scirocco_kulfan_server` runs it. Header:
[`scirocco/kulfan_discipline.h`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/include/scirocco/kulfan_discipline.h).

See [Server reference](../python/servers) to launch it.

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
convention above, non-negative coefficients and `te_upper >= te_lower` keep
the upper surface at or above the lower one. This allows zero thickness;
use thickness or area constraints to exclude degenerate designs.

## Options

Set at discipline construction.

| Option | Type | Default | Meaning |
| --- | --- | --- | --- |
| `n_cst_upper` | int | 6 | Number of upper-surface Bernstein coefficients |
| `n_cst_lower` | int | 6 | Number of lower-surface Bernstein coefficients |
| `npanel` | int (even) | 200 | Total panels; `npoint = npanel + 1` |
| `te_upper` | double | 0.0 | Fixed upper-surface TE y-offset |
| `te_lower` | double | 0.0 | Fixed lower-surface TE y-offset |

The defaults give a closed trailing edge, which is fine for inviscid
analysis. For viscous analysis set `te_upper - te_lower` to at least 1e-6
chord (the examples use 0.00126, 0.126% chord); otherwise the analysis
discipline raises an error.

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

Bounding the smoothness outputs limits the second derivative of the shape
polynomial $S(\psi)$ and discourages coefficient oscillations. This linear
constraint is not a direct bound on the curvature of the complete surface. See
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

You normally use it through `scirocco_kulfan_server` from Python rather
than constructing it in C++.

## Evaluation stations

The discipline samples the surface at cosine-clustered $\psi$ stations
on the half-chord. The number of stations is `npanel/2 + 1`, with
mirror points giving the full panel-node grid. The `psi_eval()` accessor
exposes the half-chord $\psi$ vector for testing.

## Building, launching, using

With the server running, connect from any Philote client. In Python, with
OpenMDAO as the framework:

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

- [Python overview](../python): connect the geometry and analysis servers
- [Change the shape](../python/geometry): CST coefficients from Python
- [CST design](../optimization/cst-design): inviscid and viscous optimization
