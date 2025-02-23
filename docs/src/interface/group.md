# An interface for Lie groups

```@docs
AbstractLieGroup
LieGroup
AbstractLieAlgebraTangentVector
AbstractLieGroupPoint
```

## Functions on Lie groups

```@autodocs
Modules = [LieGroups]
Pages = ["src/interface.jl"]
Order = [:function]
```

## Internal functions and macros

```@docs
LieGroups.CommonUnitarySubAlgebra
LieGroups.@default_lie_group_fallbacks
```

## A validation Lie group

```@docs
ValidationLieGroup
```

### Internal functions

```@docs
LieGroups._vLc
LieGroups._msg
```

## Literature

```@bibliography
Pages = ["group.md"]
Canonical=false
```