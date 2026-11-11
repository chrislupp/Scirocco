---
sidebar_position: 6
---

# Pseudo-transient continuation

Pure Newton on a viscous airfoil sometimes diverges. A typical pathology
is: the residual drops a few orders of magnitude, then a step lands on a
nearby root with a different transition location, the residual jumps from
$10^{-3}$ to $10^{15}$, and the next step blows up.

**Pseudo-transient continuation (PTC)** is a small modification of Newton
that can reduce this instability: it adds a damping term that increases when
the residual increases, then relaxes back toward pure Newton as the solution
converges. The target residual equations are unchanged.

PTC is **on by default** in Scirocco. This page explains what it's doing,
when it kicks in, and how to retune it if you need to.

## The modification

Standard Newton step:

$$
\mathbf{R}_U^k \,\Delta \mathbf{U} = -\mathbf{R}^k
$$

PTC step:

$$
\left(\mathbf{R}_U^k + \frac{1}{\tau^k} \mathbf{I}\right) \Delta \mathbf{U}
   = -\mathbf{R}^k
$$

Adding $1/\tau \cdot \mathbf{I}$ to the diagonal regularizes the Jacobian.
In the limit $\tau \to \infty$ it's pure Newton. In the limit $\tau \to 0$
the step approaches $-\tau R$, a residual-directed pseudo-time step.
PTC starts at large
$\tau$ (essentially Newton), shrinks $\tau$ whenever the residual jumps,
and grows $\tau$ when the residual decreases.

The converged residual is $\mathbf{R} = 0$ regardless of $\tau$, so PTC
does not change the target equations. Different paths can still reach
different roots when the nonlinear system has multiple solutions.

## The SER schedule

Scirocco updates $\tau$ each iteration with a **switched
evolution-relaxation** (SER) rule based on the residual ratio
$r^k = \|\mathbf{R}^k\| / \|\mathbf{R}^{k-1}\|$:

| Condition | Update |
| --- | --- |
| $r^k < 1$ (residual decreased) | grow $\tau$ by $\min(1/r^k,\, \text{tau\_grow\_clamp})$ |
| $r^k \ge 1$ (residual stagnant or grew) | shrink $\tau$ by $\max(1/r^k,\, \text{tau\_shrink\_clamp})$ |
| always | clamp to $[\text{tau\_min},\, \text{tau\_max}]$ |

The default constants in
[`data/param.h`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/include/scirocco/data/param.h):

| Field | Default | Role |
| --- | --- | --- |
| `use_ptc` | `true` | Master on/off |
| `tau_init` | `1e6` | Initial value — effectively Newton |
| `tau_min` | `1e-2` | Floor: maximum damping |
| `tau_max` | `1e8` | Ceiling: effectively Newton again |
| `tau_grow_clamp` | `4.0` | Per-iter growth ratio cap |
| `tau_shrink_clamp` | `0.1` | Per-iter shrink ratio cap |
| `niglob` | `200` | Newton iteration cap (raised from 50 to give PTC headroom) |

The growth/shrink clamps prevent $\tau$ from oscillating wildly in either
direction — `tau_grow_clamp = 4.0` means at most a $4\times$ increase per
step.

## When you should care

For well-behaved cases (clean airfoil, attached flow, modest $\alpha$),
PTC introduces little damping. Its convergence benefit and cost depend on
the operating point.

For rough cases — high $\alpha$, separation, transition jumping between
panels, intermediate CST iterates in an optimizer — PTC can improve
convergence, but it does not guarantee success.

### Common adjustments

| Symptom | Try |
| --- | --- |
| Newton diverges immediately at iter 1–2 | Lower `tau_init` (e.g. `1e2`); start damped from the get-go |
| Residual stalls at $10^{-6}$ | Check geometry and scaling; increase `niglob` if the residual is still falling |
| Residual oscillates between two values | Lower `tau_grow_clamp` to 2.0 to slow the relaxation back to Newton |
| You want to compare to a non-PTC baseline | `param.use_ptc = false` |

## Adjoint considerations

PTC adds $1/\tau \cdot \mathbf{I}$ to the Jacobian **only** during the
forward Newton iteration. The stored adjoint Jacobian is the unmodified
$\mathbf{R}_U$ at the converged state, regardless of the final pseudo-time step.

If both PTC-on and PTC-off runs converge to the same state, you should get
matching gradients within numerical noise.
