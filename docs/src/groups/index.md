# An overview of Lie groups

# Alphabetical list of Lie groups

| Group | Manifold | ``∘`` | Comment |
|:------|:---------|:---------:|:------|
| [`CircleGroup`](@ref) | real or complex [`Circle`](@extref `Manifolds.Circle`), [`Sphere`](@extref `Manifolds.Sphere`) | [`*`](@ref AbelianMultiplicationGroupOperation) | |
| [`GeneralLinearGroup`](@ref) | [`InvertibleMatrices`](@extref `Manifolds.InvertibleMatrices`) | [`*`](@ref MatrixMultiplicationGroupOperation) | |
| [`HeisenbergGroup`](@ref) | [`HeisenbergMatrices`](@extref `Manifolds.HeisenbergMatrices`) | [`*`](@ref MatrixMultiplicationGroupOperation) | |
| [`OrthogonalGroup`](@ref) | [`OrthogonalMatrices`](@extref `Manifolds.OrthogonalMatrices`) | [`*`](@ref MatrixMultiplicationGroupOperation) | This can be interpreted as all rotations and reflections.
| [`PowerLieGroup`](@ref) | [`PowerManifold`](@extref `ManifoldsBase.PowerManifold`) | [`∘`](@ref PowerGroupOperation) | [`^`](@ref PowerLieGroup) is a constructor |
| [`ProductLieGroup`](@ref) | [`ProductManifold`](@extref `ManifoldsBase.ProductManifold`) | [`∘`](@ref ProductGroupOperation) | [`×`](@ref LinearAlgebra.cross(::AbstractGroupOperation...)) of two Lie groups is a constructor |
| [`LeftSemidirectProductLieGroup`](@ref) | [`ProductManifold`](@extref `ManifoldsBase.ProductManifold`) | [`∘`](@ref LeftSemidirectProductGroupOperation) | [`⋉`](@ref ⋉(L1::LieGroup, L2::LieGroup)) of 2 Lie groups is a constructor, similarly [`⋊`](@ref ⋊(L1::LieGroup, L2::LieGroup)) for the right variant |
| [`SpecialEuclideanGroup`](@ref) | [`Rotations`](@extref `Manifolds.Rotations`)[`⋉`](@ref LeftSemidirectProductGroupOperation)[`Euclidean`](@extref `Manifolds.Rotations`) | [`∘`](@ref LeftSemidirectProductGroupOperation) | Analogously you can also use a [`⋊`](@ref RightSemidirectProductGroupOperation) if you prefer tuples `(t,R)` having the rotation matrix in the second component |
| [`SpecialLinearGroup`](@ref) | [`GeneralUnitaryMatrices`](@extref `Manifolds.GeneralUnitaryMatrices`) with determinant one | [`*`](@ref MatrixMultiplicationGroupOperation) | |
| [`SpecialOrthogonalGroup`](@ref) | [`Rotations`](@extref `Manifolds.Rotations`) | [`*`](@ref MatrixMultiplicationGroupOperation) | |
| [`SpecialUnitaryGroup`](@ref) | [`GeneralUnitaryMatrices`](@extref `Manifolds.GeneralUnitaryMatrices`) | [`*`](@ref MatrixMultiplicationGroupOperation) | |
| [`SymplecticGroup`](@ref) | [`SymplecticMatrices`](@extref `Manifolds.SymplecticMatrices`) | [`*`](@ref MatrixMultiplicationGroupOperation) | |
| [`TranslationGroup`](@ref) | [`Euclidean`](@extref `Manifolds.Euclidean`) | [`+`](@ref AdditionGroupOperation) | |
| [`UnitaryGroup`](@ref) | [`UnitaryMatrices`](@extref `Manifolds.UnitaryMatrices`) | [`*`](@ref MatrixMultiplicationGroupOperation) | |
