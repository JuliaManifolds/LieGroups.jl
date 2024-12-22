# The orthogonal group

```@docs
OrthogonalGroup
```

For this Lie group, several implementations are already covered by the defaults in [the generic (matrix) multiplication operation](@ref multiplication-operation-sec).

# Functions

```@autodocs
Modules = [LieGroups]
Pages = ["groups/orthogonal_group.jl"]
Order = [:function]
```

# Utility functions

```@docs
LieGroups.angles_4d_skew_sym_matrix
LieGroups.eigen_safe
LieGroups.usinc_from_cos
```