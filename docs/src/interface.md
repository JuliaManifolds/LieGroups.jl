# The interface for Lie Groups

```@docs
LieGroup
LieAlgebra
```

## Group Operations

```@docs
AbstractGroupOperation
Identity
```

## Properties

```@docs
LieGroups.AbstractInvarianceTrait
LieGroups.HasBiinvariantMetric
LieGroups.HasLeftInvariantMetric
LieGroups.HasRightInvariantMetric
```

## Functions on Lie Groups

```@docs
compose
compose!
compose_diff_left
compose_diff_left!
compose_diff_right
compose_diff_right!
compose_inv_left
compose_inv_left!
compose_inv_right
compose_inv_right!
conjugate
conjugate!
identity_element
identity_element!
is_identity
inv
inv!
inv_diff
inv_diff!
```

## Actions on Lie Groups

## Functions on Lie Algebras

```@docs
Lie_bracket
Lie_bracket!
```

## Specific Group Operations

For some generic group operations, an easy generic implementation can be provided. This section lists these specific cases and what is implemented for these.

### Additive Group operation

```@docs
AdditiveGroupOperation
```

### Multiplicative Group operation

## Literature

```@bibliography
Pages = ["interface.md"]
Canonical=false
```
