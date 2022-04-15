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
    ⋉

include("liegroup.jl")
include("rotations.jl")
include("translations.jl")

end
