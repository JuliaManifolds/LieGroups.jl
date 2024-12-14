# An overview of Lie groups

# Alphabetical list of Lie groups

| Group | Manifold | ``∘`` | Comment |
|:------|:---------|:---------:|:------|
| [`GeneralLinearGroup`](@ref) | [`InvertibleMatrices`](@extref `Manifolds.InvertibleMatrices`) | [`*`](@ref MatrixMultiplicationGroupOperation) | |
| [`SpecialEuclideanGroup`](@ref) | [`Rotations`](@extref `Manifolds.Rotations`)`×`[`Euclidean`](@extref `Manifolds.Rotations`) | [`∘`](@ref LeftSemidirectProductGroupOperation) | |
| [`SpecialOrthogonalGroup`](@ref) | [`Rotations`](@extref `Manifolds.Rotations`) | [`*`](@ref MatrixMultiplicationGroupOperation) | |
| [`TranslationGroup`](@ref) | [`Euclidean`](@extref `Manifolds.Euclidean`) | [`+`](@ref AdditionGroupOperation) | |