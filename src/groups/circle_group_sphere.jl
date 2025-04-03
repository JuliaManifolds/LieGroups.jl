#
#circle group represented as vectors in ℝ^2, operation: adding angles in the unit circle
#
function CircleGroup(M::Sphere{ManifoldsBase.TypeParameter{Tuple{1}},ℝ})
    return CircleGroup{ℝ,MatrixMultiplicationGroupOperation,typeof(M)}(
        M, MatrixMultiplicationGroupOperation()
    )
end

function CircleGroup(::Euclidean{ManifoldsBase.TypeParameter{Tuple{2}},ℝ})
    return CircleGroup(Sphere(1))
end

const _PlaneCircleGroup = LieGroup{
    ℝ,MatrixMultiplicationGroupOperation,<:Sphere{ManifoldsBase.TypeParameter{Tuple{1}},ℝ}
}
