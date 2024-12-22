# A collection of utility functions

@doc raw"""
    angles_4d_skew_sym_matrix(A)

The Lie algebra of [`Orthogonal`](@ref)`(4)` in ``ℝ^{4×4}``, ``𝔰𝔬(4)``, consists of ``4×4``
skew-symmetric matrices. The unique imaginary components of their eigenvalues are the
angles of the two plane rotations. This function computes these more efficiently than
`eigvals`.

By convention, the returned values are sorted in decreasing order.
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

@doc raw"""
    usinc_from_cos(x::Real)

Unnormalized version of `sinc` function, i.e. ``\operatorname{usinc}(θ) = \frac{\sin(θ)}{θ}``,
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
