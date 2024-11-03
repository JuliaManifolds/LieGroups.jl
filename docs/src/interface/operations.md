
# An interface for Lie group operations

```@docs
AbstractGroupOperation
Identity
```

You can combine some specific group operations with one of several manifolds to form a Lie group.
You can still define the corresponding functions generically for all groups with this group operation regardless of the manifold.
The following sections collect these.


* an [`AdditionGroupOperation`](@ref)

## [Additive group operation](@id addition-operation-sec)

```@autodocs
Modules = [LieGroups]
Pages = ["addition_operation.jl"]
Order = [:type, :function]
```

## [Multiplication group operation](@id multiplication-operation-sec)

```@autodocs
Modules = [LieGroups]
Pages = ["multiplication_operation.jl"]
Order = [:type, :function]
```

## Literature

```@bibliography
Pages = ["operations.md"]
Canonical=false
```