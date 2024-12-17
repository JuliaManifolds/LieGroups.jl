@doc raw"""
    CircleGroup <: GroupManifold{Circle{ℂ},MultiplicationOperation}

The circle group is the complex circle ([`Circle(ℂ)`](@ref)) equipped with
the group operation of complex multiplication ([`MatrixMultiplicationGroupOperation`](@ref)).
"""

const CircleGroup = LieGroup{
   ℂ, MatrixMultiplicationGroupOperation, Manifolds.AbstractSphere{ℂ}
}

function CircleGroup()
    circ = Manifolds.Sphere(0, ℂ)
    return CircleGroup(
        circ, MatrixMultiplicationGroupOperation()
        )
end