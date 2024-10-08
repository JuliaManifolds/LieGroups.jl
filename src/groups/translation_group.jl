"""
    TranslationGroup{𝔽}

The Lie group consisting of the [`AdditionGroupOperation`](@ref) on some
[`Euclidean`](@extref `Manifolds.Euclidean`) space.

# Constructor

TODO
"""
const TranslationGroup{𝔽} = LieGroup{𝔽,Euclidean,AdditionGroupOperation}
