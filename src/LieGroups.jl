module LieGroups

using LinearAlgebra

import Base: identity, +, -, *, inv, ==

export
    # liegroup
    AbstractLieGroup,
    AbstractLieAlgebra,
    dim,
    dof,

    # rotations
    AbstractRotationGroup,
    AbstractRotationAlgebra,
    SO, so,
    ∧, ∨,
    ⋉,
    jacobian,
    ⊕,

    # rigid_motions
    SpecialEuclideanGroup,
    SpecialEuclideanAlgebra,
    SE, se

include("utils.jl")
include("liegroup.jl")
include("rotations.jl")
include("rigid_motions.jl")

end
