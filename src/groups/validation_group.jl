"""
    ValidationLieGroup{L<:AbstractLieGroup} <: AbstractLieGroup

A Lie group to add tests to input parameters and ouptut values of functions defined
for [`LieGroups`](@ref).


# Fields

* `lie_group::L` the [`AbstractLieGroup`](@ref) to be decorated
* ``

# Constructor
    ValidationLieGroup(L::AbstractLieGroup; )

Generate the Validation Lie Group for the given [`AbstractLieGroup`](@ref)  `L`

    ValidationLieGroup(M::AbstractManifold, op::AbstractGroupOperation; kwargs...)

Generate the Validation Lie Group for the given [`AbstractLieGroup`](@ref)  `L` based on
a the [`ValidationManifold`](@extref) of `M` and a group operation `op`, that is, also on the
manifold operations can be tested. suitable keywords are passed down.

# Keyword arguments

For the second constructor, all further keywords are passed to the [`ValidationManifold`](@ref) as well.
"""
struct ValidationLieGroup{𝔽, L<:LieGroup{𝔽}}
    lie_group::L
    mode::S
end