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
adjoint
adjoint!
base_manifold
compose
compose!
diff_left_compose
diff_left_compose!
diff_right_compose
diff_right_compose!
inv_left_compose
inv_left_compose!
inv_right_compose
inv_right_compose!
conjugate
conjugate!
diff_conjugate
diff_conjugate!
exp
exp!
identity_element
identity_element!
is_identity
inv(::LieGroup, ::Any)
inv!
diff_inv
diff_inv!
log
log!
```

## Actions on Lie groups

```@autodocs
Modules = [LieGroups]
Pages = ["group_action_interface.jl"]
Order = [:type]
```

### Functions for Lie group actions

```@autodocs
Modules = [LieGroups]
Pages = ["group_action_interface.jl"]
Order = [:function]
```

### Specific Lie group actions

```@autodocs
Modules = [LieGroups]
Pages = ["group_operation_action.jl"]
Order = [:type, :function]
```

## Functions on Lie Algebras

```@docs
lie_bracket
lie_bracket!
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
