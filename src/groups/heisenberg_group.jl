
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

function _heisenberg_a_view(G::HeisenbergGroup, g)
    n = ManifoldsBase.get_parameter(G.manifold.size)[1]
    return view(g, 1, 2:(n + 1))
end
function _heisenberg_b_view(G::HeisenbergGroup, g)
    n = ManifoldsBase.get_parameter(G.manifold.size)[1]
    return view(g, 2:(n + 1), n + 2)
end

@doc raw"""
    exp(G::HeisenbergGroup, ::Identity{MatrixMultiplicationGroupOperation}, X)

Lie group exponential for the [`HeisenbergGroup`](@ref) `G` of the vector `X`.
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
function Base.exp(G::HeisenbergGroup, X)
    h = similar(X)
    exp!(G, h, X)
    return h
end

function ManifoldsBase.exp!(G::HeisenbergGroup, h, X)
    n = ManifoldsBase.get_parameter(G.manifold.size)[1]
    copyto!(h, I)
    a_view = _heisenberg_a_view(G, X)
    b_view = _heisenberg_b_view(G, X)
    h[1, 2:(n + 1)] .= a_view
    h[2:(n + 1), n + 2] .= b_view
    h[1, n + 2] = X[1, n + 2] + dot(a_view, b_view) / 2
    return h
end

@doc raw"""
    exp(G::HeisenbergGroup, g, X)

Exponential map on the [`HeisenbergGroup`](@ref) `G` with the left-invariant metric.
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
function Base.exp(G::HeisenbergGroup, g, X)
    h = similar(X)
    exp!(G, h, g, X)
    return h
end

function ManifoldsBase.exp!(G::HeisenbergGroup, h, g, X)
    n = ManifoldsBase.get_parameter(G.manifold.size)[1]
    copyto!(h, I)
    a_p_view = _heisenberg_a_view(G, g)
    b_p_view = _heisenberg_b_view(G, g)
    a_X_view = _heisenberg_a_view(G, X)
    b_X_view = _heisenberg_b_view(G, X)
    h[1, 2:(n + 1)] .= a_p_view .+ a_X_view
    h[2:(n + 1), n + 2] .= b_p_view .+ b_X_view
    h[1, n + 2] =
        g[1, n + 2] + X[1, n + 2] + dot(a_X_view, b_X_view) / 2 + dot(a_p_view, b_X_view)
    return h
end

@doc raw"""
    injectivity_radius(G::HeisenbergGroup)

Return the injectivity radius on the [`HeisenbergGroup`](@ref) `G`, which is ``∞``.
"""
ManifoldsBase.injectivity_radius(::HeisenbergGroup) = Inf

@doc raw"""
    log(G::HeisenbergGroup, g, h)

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
Base.log(::HeisenbergGroup, g, h)

function ManifoldsBase.log!(G::HeisenbergGroup, X, g, h)
    n = ManifoldsBase.get_parameter(G.manifold.size)[1]
    fill!(X, 0)
    a_p_view = _heisenberg_a_view(G, g)
    b_p_view = _heisenberg_b_view(G, g)
    a_q_view = _heisenberg_a_view(G, h)
    b_q_view = _heisenberg_b_view(G, h)
    X[1, 2:(n + 1)] .= a_q_view .- a_p_view
    X[2:(n + 1), n + 2] .= b_q_view .- b_p_view
    pinvq_c = dot(a_p_view, b_p_view) - g[1, n + 2] + h[1, n + 2] - dot(a_p_view, b_q_view)
    X[1, n + 2] = pinvq_c - dot(a_q_view - a_p_view, b_q_view - b_p_view) / 2
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

@doc raw"""
    log(G::HeisenbergGroup, g)

Lie group logarithm for the [`HeisenbergGroup`](@ref) `G` of the point `g`.
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
log(G::HeisenbergGroup, g)

function ManifoldsBase.log!(G::HeisenbergGroup, X, g)
    n = ManifoldsBase.get_parameter(G.manifold.size)[1]
    fill!(X, 0)
    view_a_X = _heisenberg_a_view(G, X)
    view_b_X = _heisenberg_b_view(G, X)
    view_a_X .= _heisenberg_a_view(G, g)
    view_b_X .= _heisenberg_b_view(G, g)
    X[1, n + 2] = g[1, n + 2] - dot(view_a_X, view_b_X) / 2
    return X
end
function ManifoldsBase.log!(
    ::HeisenbergGroup, X, ::Identity{MatrixMultiplicationGroupOperation}
)
    fill!(X, 0)
    return X
end

function Base.show(
    io::IO, ::HeisenbergGroup{ManifoldsBase.TypeParameter{Tuple{n}}}
) where {n}
    return print(io, "HeisenbergGroup($(n))")
end
function Base.show(io::IO, G::HeisenbergGroup{Tuple{Int}})
    n = ManifoldsBase.get_parameter(G.manifold.size)[1]
    return print(io, "HeisenbergGroup($(n); parameter=:field)")
end
