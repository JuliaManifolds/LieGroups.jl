
@doc raw"""
    HeisenbergGroup{T} <: AbstractDecoratorManifold{ℝ}

Heisenberg group `HeisenbergGroup(n)` is the group of ``(n+2)×(n+2)`` matrices [BinzPods:2008](@cite)

```math
\begin{bmatrix} 1 & \mathbf{a} & c \\
\mathbf{0} & I_n & \mathbf{b} \\
0 & \mathbf{0} & 1 \end{bmatrix}
```

where ``I_n`` is the ``n×n`` unit matrix, ``\mathbf{a}`` is a row vector of length ``n``,
``\mathbf{b}`` is a column vector of length ``n`` and ``c`` is a real number.
The group operation is matrix multiplication.

The left-invariant metric on the manifold is used.
"""
const HeisenbergGroup{T} = LieGroup{
    ℝ,MatrixMultiplicationGroupOperation,Manifolds.HeisenbergMatrices{T}
}

function HeisenbergGroup(n::Int; parameter::Symbol=:type)
    Hm = Manifolds.HeisenbergMatrices(n; parameter=parameter)
    return HeisenbergGroup{typeof(Hm).parameters...}(
        Hm, MatrixMultiplicationGroupOperation()
    )
end

function _heisenberg_a_view(M::HeisenbergGroup, p)
    n = ManifoldsBase.get_parameter(M.manifold.size)[1]
    return view(p, 1, 2:(n + 1))
end
function _heisenberg_b_view(M::HeisenbergGroup, p)
    n = ManifoldsBase.get_parameter(M.manifold.size)[1]
    return view(p, 2:(n + 1), n + 2)
end

embed(::HeisenbergGroup, p) = p
embed(::HeisenbergGroup, p, X) = X

@doc raw"""
    exp(M::HeisenbergGroup, ::Identity{MatrixMultiplicationGroupOperation}, X)

Lie group exponential for the [`HeisenbergGroup`](@ref) `M` of the vector `X`.
The formula reads
```math
\exp\left(\begin{bmatrix} 0 & \mathbf{a} & c \\
\mathbf{0} & 0_n & \mathbf{b} \\
0 & \mathbf{0} & 0 \end{bmatrix}\right) = \begin{bmatrix} 1 & \mathbf{a} & c + \mathbf{a}⋅\mathbf{b}/2 \\
\mathbf{0} & I_n & \mathbf{b} \\
0 & \mathbf{0} & 1 \end{bmatrix}
```
where ``I_n`` is the ``n×n`` identity matrix, ``0_n`` is the ``n×n`` zero matrix
and ``\mathbf{a}⋅\mathbf{b}`` is dot product of vectors.
"""
function Base.exp(M::HeisenbergGroup, e::Identity{MatrixMultiplicationGroupOperation}, X)
    q = similar(X)
    exp!(M, q, e, X)
    return q
end

function ManifoldsBase.exp!(
    M::HeisenbergGroup, q, ::Identity{MatrixMultiplicationGroupOperation}, X
)
    n = ManifoldsBase.get_parameter(M.manifold.size)[1]
    copyto!(q, I)
    a_view = _heisenberg_a_view(M, X)
    b_view = _heisenberg_b_view(M, X)
    q[1, 2:(n + 1)] .= a_view
    q[2:(n + 1), n + 2] .= b_view
    q[1, n + 2] = X[1, n + 2] + dot(a_view, b_view) / 2
    return q
end

@doc raw"""
    exp(M::HeisenbergGroup, p, X)

Exponential map on the [`HeisenbergGroup`](@ref) `M` with the left-invariant metric.
The expression reads
```math
\exp_{\begin{bmatrix} 1 & \mathbf{a}_p & c_p \\
\mathbf{0} & I_n & \mathbf{b}_p \\
0 & \mathbf{0} & 1 \end{bmatrix}}\left(\begin{bmatrix} 0 & \mathbf{a}_X & c_X \\
\mathbf{0} & 0_n & \mathbf{b}_X \\
0 & \mathbf{0} & 0 \end{bmatrix}\right) =
\begin{bmatrix} 1 & \mathbf{a}_p + \mathbf{a}_X & c_p + c_X + \mathbf{a}_X⋅\mathbf{b}_X/2 + \mathbf{a}_p⋅\mathbf{b}_X \\
\mathbf{0} & I_n & \mathbf{b}_p + \mathbf{b}_X \\
0 & \mathbf{0} & 1 \end{bmatrix}
```
where ``I_n`` is the ``n×n`` identity matrix, ``0_n`` is the ``n×n`` zero matrix
and ``\mathbf{a}⋅\mathbf{b}`` is dot product of vectors.
"""
function Base.exp(M::HeisenbergGroup, p, X)
    q = similar(X)
    exp!(M, q, p, X)
    return q
end

function ManifoldsBase.exp!(M::HeisenbergGroup, q, p, X)
    n = ManifoldsBase.get_parameter(M.manifold.size)[1]
    copyto!(q, I)
    a_p_view = _heisenberg_a_view(M, p)
    b_p_view = _heisenberg_b_view(M, p)
    a_X_view = _heisenberg_a_view(M, X)
    b_X_view = _heisenberg_b_view(M, X)
    q[1, 2:(n + 1)] .= a_p_view .+ a_X_view
    q[2:(n + 1), n + 2] .= b_p_view .+ b_X_view
    q[1, n + 2] =
        p[1, n + 2] + X[1, n + 2] + dot(a_X_view, b_X_view) / 2 + dot(a_p_view, b_X_view)
    return q
end

@doc raw"""
    injectivity_radius(M::HeisenbergGroup)

Return the injectivity radius on the [`HeisenbergGroup`](@ref) `M`, which is ``∞``.
"""
injectivity_radius(::HeisenbergGroup) = Inf

@doc raw"""
    log(G::HeisenbergGroup, p, q)

Compute the logarithmic map on the [`HeisenbergGroup`](@ref) group.
The formula reads
```math
\log_{\begin{bmatrix} 1 & \mathbf{a}_p & c_p \\
\mathbf{0} & I_n & \mathbf{b}_p \\
0 & \mathbf{0} & 1 \end{bmatrix}}\left(\begin{bmatrix} 1 & \mathbf{a}_q & c_q \\
\mathbf{0} & I_n & \mathbf{b}_q \\
0 & \mathbf{0} & 1 \end{bmatrix}\right) =
\begin{bmatrix} 0 & \mathbf{a}_q - \mathbf{a}_p & c_q - c_p + \mathbf{a}_p⋅\mathbf{b}_p - \mathbf{a}_q⋅\mathbf{b}_q - (\mathbf{a}_q - \mathbf{a}_p)⋅(\mathbf{b}_q - \mathbf{b}_p) / 2 \\
\mathbf{0} & 0_n & \mathbf{b}_q - \mathbf{b}_p \\
0 & \mathbf{0} & 0 \end{bmatrix}
```
where ``I_n`` is the ``n×n`` identity matrix, ``0_n`` is the ``n×n`` zero matrix
and ``\mathbf{a}⋅\mathbf{b}`` is dot product of vectors.
"""
Base.log(::HeisenbergGroup, p, q)

function ManifoldsBase.log!(M::HeisenbergGroup, X, p, q)
    n = ManifoldsBase.get_parameter(M.manifold.size)[1]
    fill!(X, 0)
    a_p_view = _heisenberg_a_view(M, p)
    b_p_view = _heisenberg_b_view(M, p)
    a_q_view = _heisenberg_a_view(M, q)
    b_q_view = _heisenberg_b_view(M, q)
    X[1, 2:(n + 1)] .= a_q_view .- a_p_view
    X[2:(n + 1), n + 2] .= b_q_view .- b_p_view
    pinvq_c = dot(a_p_view, b_p_view) - p[1, n + 2] + q[1, n + 2] - dot(a_p_view, b_q_view)
    X[1, n + 2] = pinvq_c - dot(a_q_view - a_p_view, b_q_view - b_p_view) / 2
    return X
end

@doc raw"""
    log(M::HeisenbergGroup, ::Identity{MatrixMultiplicationGroupOperation}, p)

Lie group logarithm for the [`HeisenbergGroup`](@ref) `M` of the point `p`.
The formula reads
```math
\log\left(\begin{bmatrix} 1 & \mathbf{a} & c \\
\mathbf{0} & I_n & \mathbf{b} \\
0 & \mathbf{0} & 1 \end{bmatrix}\right) =
\begin{bmatrix} 0 & \mathbf{a} & c - \mathbf{a}⋅\mathbf{b}/2 \\
\mathbf{0} & 0_n & \mathbf{b} \\
0 & \mathbf{0} & 0 \end{bmatrix}
```
where ``I_n`` is the ``n×n`` identity matrix, ``0_n`` is the ``n×n`` zero matrix
and ``\mathbf{a}⋅\mathbf{b}`` is dot product of vectors.
"""
log(M::HeisenbergGroup, ::Identity{MatrixMultiplicationGroupOperation}, p)

function ManifoldsBase.log!(
    M::HeisenbergGroup, X, ::Identity{MatrixMultiplicationGroupOperation}, p
)
    n = ManifoldsBase.get_parameter(M.manifold.size)[1]
    fill!(X, 0)
    view_a_X = _heisenberg_a_view(M, X)
    view_b_X = _heisenberg_b_view(M, X)
    view_a_X .= _heisenberg_a_view(M, p)
    view_b_X .= _heisenberg_b_view(M, p)
    X[1, n + 2] = p[1, n + 2] - dot(view_a_X, view_b_X) / 2
    return X
end
function ManifoldsBase.log!(
    ::HeisenbergGroup,
    X,
    ::Identity{MatrixMultiplicationGroupOperation},
    ::Identity{MatrixMultiplicationGroupOperation},
)
    fill!(X, 0)
    return X
end

function Base.show(
    io::IO, ::HeisenbergGroup{ManifoldsBase.TypeParameter{Tuple{n}}}
) where {n}
    return print(io, "HeisenbergGroup($(n))")
end
function Base.show(io::IO, M::HeisenbergGroup{Tuple{Int}})
    n = ManifoldsBase.get_parameter(M.manifold.size)[1]
    return print(io, "HeisenbergGroup($(n); parameter=:field)")
end
