
# An interface for Lie group operations

```@docs
AbstractGroupOperation
Identity
```

Some specific group operations allow to define functions for all manifolds they
do build a Lie group with. These are

* an [`AdditionGroupOperation`](@ref)

## Additive group operation

```@autodocs
Modules = [LieGroups]
Pages = ["addition.jl"]
Order = [:type, :function]
```