module LieGroup

using LinearAlgebra

import Base: identity, +, -, *, inv, ==

export
    # liegroup
    AbstractLieGroup,
    AbstractLieAlgebra,

    # rotations
    AbstractRotationGroup,
    AbstractRotationAlgebra,
    SO, so,
    ∧, ∨,
    ⋉,

    # rigid_motions
    SpecialEuclideanGroup,
    SpecialEuclideanAlgebra,
    SE, se

include("utils.jl")
include("liegroup.jl")
include("rotations.jl")
include("rigid_motions.jl")

end
