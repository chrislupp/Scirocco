---
sidebar_position: 1
---

# Theory

Mathematical and derivation notes that go below the level of the
[Concepts](../concepts/panel-method) section.

## Adjoint sensitivity

Derivation of the discrete adjoint for the coupled inviscid + viscous
system, including the AIC structure, the dependency chain through
$\delta^*$, and the cost-counting argument for "all gradients in one
adjoint solve":

→ [`theory/adjoint.md`](https://github.com/chrislupp/Scirocco/blob/develop/theory/adjoint.md)
on GitHub.

The document lives in the repository rather than on this site so it can
stay close to the source of truth (it's edited by the same hands that
write the adjoint code, and it references file/line numbers in
`scirocco/src/`). If you find an inconsistency between the theory note and
the implementation, please file an issue against the repo.

## Coming later

The current theory writeup is focused on adjoints. Possible future
additions, not yet written:

- Discretization details of the panel method (linear vortex influence integrals)
- Integral boundary-layer derivations (momentum, shape-parameter, $C_\tau$)
- Wake closure and Squire–Young extrapolation
- $e^N$ envelope and Drela's closure fits

If you want one of these in particular, an issue on the repo helps
prioritize.
