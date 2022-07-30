module LieGroups

using LinearAlgebra

import Base: identity, +, -, *, inv, ==

export
    # liegroup
    AbstractLieGroup,
    AbstractLieAlgebra,
    dim,
    dof,
    ∧, ∨,
    ⋉,
    jacobian,
    ⊕,

    # rotations
    AbstractRotationGroup,
    AbstractRotationAlgebra,
    SO, so,
    rotation,

    # rigid_motions
    SpecialEuclideanGroup,
    SpecialEuclideanAlgebra,
    SE, se,
    translation

include("utils.jl")
include("liegroup.jl")
include("liealgebra.jl")
include("rotations.jl")
include("rigid_motions.jl")

end
