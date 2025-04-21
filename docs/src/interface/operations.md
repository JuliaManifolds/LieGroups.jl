
# An interface for Lie group operations

```@docs
AbstractGroupOperation
Identity
```

You can combine some specific group operations with one of several manifolds to form a Lie group.
You can still define the corresponding functions generically for all groups with this group operation regardless of the manifold.
The following sections collect these.


* an [`AdditionGroupOperation`](@ref)
* a [`AbelianMultiplicationGroupOperation`](@ref) and
* a [`MatrixMultiplicationGroupOperation`](@ref)


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

## [Abelian multiplication group operation](@id multiplication-operationabelian-sec)
Since the commutativity of some Lie groups leads to a significant simplification of computations on those groups, the [`abelian multiplication group operation`](@ref multiplication-operationabelian-sec) optimizes these cases.
Additionally, it provides an interface for the abelian Lie groups. Some of these can are represented by `isbits`-types,
which don't have mutating variants of the functions.

```@autodocs
Modules = [LieGroups]
Pages = ["multiplication_operation_abelian.jl"]
Order = [:type, :function]
```

## Literature

```@bibliography
Pages = ["operations.md"]
Canonical=false
```