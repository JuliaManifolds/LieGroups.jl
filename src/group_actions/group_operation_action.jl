"""
    LeftGroupOperationAction <: AbstractLeftGroupActionType

A type for the [`AbstractLeftGroupActionType`](@ref) when acting on the group itself
from the left, that is

```math
σ_h(g) = σ(h,g) = h$(_math(:∘))g
```

for its inverse ``(σ_h)^{-1}`` see [`InverseLeftGroupOperationActionAction`](@ref).
"""
struct LeftGroupOperationAction <: AbstractLeftGroupActionType end

"""
    RightGroupOperationAction <: AbstractRightGroupActionType

A type for the [`AbstractLeftGroupActionType`](@ref) when acting on the group itself
gfrom the right

```math
σ_h(g) = σ(h,g) = g$(_math(:∘))h
```

for its inverse ``(σ_h)^{-1}`` see [`InverseRightGroupOperationActionAction`](@ref).
"""
struct RightGroupOperationAction <: AbstractLeftGroupActionType end

"""
    InverseLeftGroupOperationActionAction <: AbstractRightGroupActionType

A type for the [`AbstractLeftGroupActionType`](@ref) when acting on the group itself
given by the inverse of a [`LeftGroupOperationAction`](@ref) ``σ_h`` as

```math
τ_h(g) $(_tex(:def)) σ_h^{-1}(g) = σ(h^{-1},g) = h^{-1}$(_math(:∘))g
```

Note that while in naming it is the inverse of the left action, it's
properties yield that is is an [`AbstractRightGroupActionType`](@ref), since

$(_note(:LeftInverseActionIsRight))

for its inverse ``(σ_h)^{-1}`` see [`InverseLeftGroupOperationActionAction`](@ref).

!!! note
    Some literature also calls this by itself _the_ right group operation action.
"""
struct InverseLeftGroupOperationActionAction <: AbstractRightGroupActionType end

"""
    InverseRightGroupOperationActionAction <: AbstractLeftGroupActionType

A type for the [`AbstractLeftGroupActionType`](@ref) when acting on the group itself
given by the inverse of a [`RightGroupOperationAction`](@ref) ``σ_h`` as

```math
τ_h(g) $(_tex(:def)) σ_h^{-1}(g) = σ(h^{-1},g) = g$(_math(:∘))h^{-1}
```

Note that while in naming it is the inverse of the right action, it's
properties yield that is is an [`AbstractLeftGroupActionType`](@ref), since

$(_note(:RightInverseActionIsLeft))

for its inverse ``(σ_h)^{-1}`` see [`InverseLeftGroupOperationActionAction`](@ref).
"""
struct InverseRightGroupOperationActionAction <: AbstractLeftGroupActionType end


_doc_apply_groupop = """
    apply( (G,A), g, h)

    apply(A::GroupOperationAction, g, h)
    apply!(A::GroupOperationAction, k, g, h)

apply the stored group operation action, using [`compose`](@ref) ``$(_math(:∘))``.
this can be done in-place of `k`.
"""

@doc "$(_doc_apply_groupop)"
apply(A::GroupOperationAction, g, h)

@doc "$(_doc_apply_groupop)"
apply!(A::GroupOperationAction, k, g, h)

function apply!(A::GroupOperationAction{LeftGroupOperationAction}, k, g, h)
    return compose!(A.group, k, g, h) #apply/compose g from left
end
function apply!(A::GroupOperationAction{RightGroupOperationAction}, k, g, h)
    return compose!(A.group, k, h, g) #apply/compose g from right
end
function apply!(A::GroupOperationAction{InverseLeftGroupOperationActionAction}, k, g, h)
    return inv_left_compose!(A.group, k, g, h) #apply/compose inv(g) from left
end
function apply!(A::GroupOperationAction{InverseRightGroupOperationActionAction}, k, g, h)
    return inv_right_compose!(A.group, k, g, h) #apply/compose inv(g) from right
end

_doc_diff_apply_groupop = """

    diff_apply(A::GroupOperationAction, g, p, X)

For the group operation action ``σ_g(p)``, compute the differential
``D_p σ_g(p): T_p$(_math(:G)) → T_{σ_g(p)}$(_math(:G))``, that is

* for the [`LeftGroupOperationAction`](@ref) this calls [`diff_right_compose`](@ref)`(G, g, p, X)`, since here ``σ_g(p) = g$(_math(:∘))p``
* for the [`RightGroupOperationAction`](@ref) this calls [`diff_left_compose`](@ref)`(G, p, g, X)`, since here ``σ_g(p) = p$(_math(:∘))g``
* for the [`InverseLeftGroupOperationActionAction`](@ref) this calls [`diff_right_compose`](@ref) with ``g^{-1}``, since here ``σ_g(p) = g^{-1}$(_math(:∘))p``
* for the [`InverseRightGroupOperationActionAction`](@ref) this calls [`diff_left_compose`](@ref) with ``g^{-1}``, since here ``σ_g(p) = p$(_math(:∘))g^{-1}``
"""
@doc "$(_doc_diff_apply_groupop)"
diff_apply(A::GroupOperationAction, g, p, X)

@doc "$(_doc_diff_apply_groupop)"
diff_apply!(A::GroupOperationAction, Y, g, p, X)

function diff_apply!(A::GroupOperationAction{LeftGroupOperationAction}, Y, g, p, X)
    return diff_right_compose!(A.group, Y, g, p, X)
end
function diff_apply!(A::GroupOperationAction{RightGroupOperationAction}, Y, g, p, X)
    return diff_left_compose!(A.group, Y, p, g, X)
end
function diff_apply!(A::GroupOperationAction{InverseLeftGroupOperationActionAction}, Y, g, p, X)
    return diff_right_compose!(A.group, Y, inv(A.group, g), p, X)
end
function diff_apply!(A::GroupOperationAction{InverseRightGroupOperationActionAction}, Y, g, p, X)
    return diff_left_compose!(A.group, Y, p, inv(A.group, g), X)
end

_doc_diff_group_apply_groupop = """

    diff_group_apply(A::GroupOperationAction, g, p, X)

Compute the differential ``D_g σ_g(p): $(_math(:𝔤)) → $(_math(:𝔤))`` of a group operation action,
that is

* for the [`LeftGroupOperationAction`](@ref) this calls [`diff_left_compose`](@ref)`(G, g, p, X)`, since here ``σ_g(p) = g$(_math(:∘))p``
* for the [`RightGroupOperationAction`](@ref) this calls [`diff_right_compose`](@ref)`(G, p, g, X)`, since here ``σ_g(p) = p$(_math(:∘))g``
* for the [`InverseLeftGroupOperationActionAction`](@ref) this calls [`diff_left_compose`](@ref) with ``g^{-1}``, since here ``σ_g(p) = g^{-1}$(_math(:∘))p`` together with [`diff_inv`](@ref)
* for the [`InverseRightGroupOperationActionAction`](@ref) this calls [`diff_right_compose`](@ref) with ``g^{-1}``, since here ``σ_g(p) = p$(_math(:∘))g^{-1}`` together with [`diff_inv`](@ref)
"""

@doc "$(_doc_diff_apply_groupop)"
diff_group_apply(A::GroupOperationAction, g, p, X)

@doc "$(_doc_diff_apply_groupop)"
diff_group_apply!(A::GroupOperationAction, Y, g, p, X)

function diff_group_apply!(A::GroupOperationAction{LeftGroupOperationAction}, Y, g, p, X)
    return diff_left_compose!(A.group, Y, g, p, X)
end
function diff_group_apply!(A::GroupOperationAction{RightGroupOperationAction}, Y, g, p, X)
    return diff_right_compose!(A.group, Y, p, g, X)
end
function diff_group_apply!(A::GroupOperationAction{InverseLeftGroupOperationActionAction}, Y, g, p, X)
    diff_inv!(A.group, Y, g, X)
    return diff_left_compose!(A.group, Y, inv(A.group, g), p, Y)
end
function diff_group_apply!(A::GroupOperationAction{InverseRightGroupOperationActionAction}, Y, g, p, X)
    diff_inv!(A.group, Y, g, X)
    return diff_right_compose!(A.group, Y, p, inv(A.group, g), Y)
end

"""
    inv(::GroupOperationAction)

return the inverse group operation action, that is, use the type representing the
inverse operation.
"""
Base.inv(A::GroupOperationAction) = GroupOperationAction(inv(A.type), A.group)

"""
    inv(::LeftGroupOperationAction)

Return the inverse of the [`LeftGroupOperationAction`](@ref), that is the [`InverseLeftGroupOperationActionAction`](@ref).
"""
Base.inv(::LeftGroupOperationAction) = InverseLeftGroupOperationActionAction()
"""
    inv(::RightGroupOperationAction)

Return the inverse of the [`RightGroupOperationAction`](@ref), that is the [`InverseRightGroupOperationActionAction`](@ref).
"""
Base.inv(::RightGroupOperationAction) = InverseRightGroupOperationActionAction()
"""
    inv(::InverseLeftGroupOperationActionAction)

Return the inverse of the [`InverseLeftGroupOperationActionAction`](@ref), that is the [`LeftGroupOperationAction`](@ref).
"""
Base.inv(::InverseLeftGroupOperationActionAction) = LeftGroupOperationAction()
"""
    inv(::InverseRightGroupOperationActionAction)

Return the inverse of the [`InverseRightGroupOperationActionAction`](@ref), that is the [`RightGroupOperationAction`](@ref).
"""
Base.inv(::InverseRightGroupOperationActionAction) = RightGroupOperationAction()

function Base.show(io::IO, A::GroupOperationAction)
    return print(io, "GroupOperationAction($(A.type), $(A.group))")
end

"""
    switch(::LeftGroupOperationAction)

Return the [`RightGroupOperationAction`](@ref), that is,
turns ``σ_g = g$(_math(:∘))h`` into ``τ_g(h) = h$(_math(:∘))g``
"""
switch(::LeftGroupOperationAction) = RightGroupOperationAction()

"""
    switch(::RightGroupOperationAction)

Return the [`LeftGroupOperationAction`](@ref), that is,
turns ``σ_g = h$(_math(:∘))g`` into ``τ_g(h) = g$(_math(:∘))h``
"""
switch(::RightGroupOperationAction) = LeftGroupOperationAction()

"""
    switch(::InverseLeftGroupOperationActionAction)

Return the [`InverseRightGroupOperationActionAction`](@ref), that is,
turns ``σ_g = g^{-1}$(_math(:∘))h`` into ``τ_g(h) = h$(_math(:∘))g^{-1}``
"""
switch(::InverseLeftGroupOperationActionAction) = InverseRightGroupOperationActionAction()

"""
    switch(::InverseRightGroupOperationActionAction)

Return the [`InverseLeftGroupOperationActionAction`](@ref), that is,
turns ``σ_g = h$(_math(:∘))g^{-1}`` into ``τ_g(h) = g^{-1}$(_math(:∘))h``
"""
switch(::InverseRightGroupOperationActionAction) = InverseLeftGroupOperationActionAction()
