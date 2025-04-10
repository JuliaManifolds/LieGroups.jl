#
#circle group represented as vectors in ℝ^2, operation: adding angles in the unit circle
#
function CircleGroup(M::Sphere{ManifoldsBase.TypeParameter{Tuple{1}},ℝ})
    return LieGroup{ℝ,AbelianMultiplicationGroupOperation,typeof(M)}(
        M, AbelianMultiplicationGroupOperation()
    )
end

function CircleGroup(::Euclidean{ManifoldsBase.TypeParameter{Tuple{2}},ℝ})
    return CircleGroup(Sphere(1))
end

const _PlanarCircleGroup = LieGroup{
    ℝ,AbelianMultiplicationGroupOperation,<:Sphere{ManifoldsBase.TypeParameter{Tuple{1}},ℝ}
}

identity_element(G::_PlanarCircleGroup) = [1.0,0.0]
function identity_element(::_PlanarCircleGroup, p::Array{<:Any, 1})
    p = [1.0,0.0]
    return p
end
function _compose(::_PlanarCircleGroup, p, q)
    a = p[1]
    b = p[2]
    c = q[1]
    d = q[2]
    z = [(a*c - b*d),(a*d + b*c)]
    return z
end
Base.inv(::_PlanarCircleGroup, p) = [p[1], -p[2]]

function inv_left_compose(C::_PlanarCircleGroup, g, h)
    g1 = Base.inv(C, g)
    return compose(C, g1, h)
end

function inv_right_compose(C::_PlanarCircleGroup, g, h)
    h1 = Base.inv(C, h)
    return compose(C, g, h1)
end

function Base.show(io::IO, ::_PlanarCircleGroup)
    return print(io, "CircleGroup(Sphere(1))")
end