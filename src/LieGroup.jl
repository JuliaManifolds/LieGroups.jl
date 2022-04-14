module LieGroup

import Base: identity, +, -, *, inv, ==

export 
    # rotations
    AbstractRotationGroup,
    IdentityRotationGroup,
    SO

include("rotations.jl")
include("translations.jl")

end
