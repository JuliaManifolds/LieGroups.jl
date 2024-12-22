# A collection of utility functions

@doc raw"""
    angles_4d_skew_sym_matrix(A)

The Lie algebra ``𝔰𝔬(4)`` of [`OrthogonalGroup`](@ref)`(4)` in ``ℝ^{4×4}``, consists of ``4×4``
skew-symmetric matrices. The unique imaginary components of their eigenvalues are the
angles of the two plane rotations. This function computes these more efficiently than
`eigvals`.

By convention, the returned values are sorted in decreasing order.
See also [`cos_angles_4d_rotation_matrix`](@ref).
"""
function angles_4d_skew_sym_matrix(A)
    @assert size(A) == (4, 4)
    @inbounds begin
        halfb = (A[1, 2]^2 + A[1, 3]^2 + A[2, 3]^2 + A[1, 4]^2 + A[2, 4]^2 + A[3, 4]^2) / 2
        c = (A[1, 2] * A[3, 4] - A[1, 3] * A[2, 4] + A[1, 4] * A[2, 3])^2
    end
    sqrtdisc = sqrt(halfb^2 - c)
    return sqrt(halfb + sqrtdisc), sqrt(halfb - sqrtdisc)
end

_doc_cos_amlges_4d_rotation_matrix = """
    cos_angles_4d_rotation_matrix(R)

4D rotations can be described by two orthogonal planes that are unchanged by
the action of the rotation (vectors within a plane rotate only within the
plane). The cosines of the two angles ``α,β`` of rotation about these planes may be
obtained from the distinct real parts of the eigenvalues of the rotation
matrix. This function computes these more efficiently by solving the system

```math
$(_tex(:aligned,
    "$(_tex(:cos))α + $(_tex(:cos))β &= $(_tex(:frac, "1","2"))$(_tex(:rm,"tr"))(R)",
    "$(_tex(:cos))α $(_tex(:cos))β &= $(_tex(:frac, "1","8"))$(_tex(:rm,"tr"))(R)^2 - $(_tex(:frac, "1","16"))$(_tex(:rm,"tr"))((R - R^T)^2) - 1"
))
```

By convention, the returned values are sorted in decreasing order.
See also [`angles_4d_skew_sym_matrix`](@ref).
"""

@doc "$(_doc_cos_amlges_4d_rotation_matrix)"
function cos_angles_4d_rotation_matrix(R)
    a = tr(R)
    b = sqrt(clamp(2 * dot(transpose(R), R) - a^2 + 8, 0, Inf))
    return ((a + b) / 4, (a - b) / 4)
end

"""
    eigen_safe(x)

Compute the eigendecomposition of `x`. If `x` is a `StaticMatrix`, it is
converted to a `Matrix` before the decomposition.
"""
@inline eigen_safe(x; kwargs...) = eigen(x; kwargs...)
@inline function eigen_safe(x::StaticMatrix; kwargs...)
    s = size(x)
    E = eigen!(Matrix(parent(x)); kwargs...)
    return Eigen(SizedVector{s[1]}(E.values), SizedMatrix{s...}(E.vectors))
end

"""
    log_safe(x)

Compute the matrix logarithm of `x`. If `x` is a `StaticMatrix`, it is
converted to a `Matrix` before computing the log.
"""
@inline log_safe(x) = log(x)
@inline function log_safe(x::StaticMatrix)
    s = Size(x)
    return SizedMatrix{s[1],s[2]}(log(Matrix(parent(x))))
end

@doc """
    usinc_from_cos(x::Real)

Unnormalized version of `sinc` function, i.e. ``$(_tex(:rm, "usinc"))(θ) = $(_tex(:frac, "$(_tex(:sin))(θ)", "θ"))``,
computed from ``x = cos(θ)``.
"""
@inline function usinc_from_cos(x::Real)
    return if x >= 1
        one(x)
    elseif x <= -1
        zero(x)
    else
        sqrt(1 - x^2) / acos(x)
    end
end
