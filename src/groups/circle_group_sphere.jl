#
#circle group represented as vectors in ℝ^2, operation: adding angles in the unit circle
#
function CircleGroup(M::Sphere{ManifoldsBase.TypeParameter{Tuple{1}},ℝ})
    return LieGroup{ℝ,AbelianMultiplicationGroupOperation,typeof(M)}(
        M,
        AbelianMultiplicationGroupOperation(),
    )
end

function CircleGroup(::Euclidean{ManifoldsBase.TypeParameter{Tuple{2}},ℝ})
    return CircleGroup(Sphere(1))
end

const _PlaneCircleGroup = LieGroup{
    ℝ,
    AbelianMultiplicationGroupOperation,
    <:Sphere{ManifoldsBase.TypeParameter{Tuple{1}},ℝ},
}
