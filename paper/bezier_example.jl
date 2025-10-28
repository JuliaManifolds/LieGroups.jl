using Pkg; Pkg.activate(@__DIR__); Pkg.instantiate()
using LieGroups, ManifoldsBase, RecursiveArrayTools
using CairoMakie
"""
    deCasteljauLieGroup(G, cp, t)

Evaluate the Beziér curve at parameter `t` given
a vector of control points `cp` on the Lie group `G`.
"""
function deCasteljauLieGroup(
        G::AbstractLieGroup, cp::Vector{P}, t::Real
    ) where {P}
    n = length(cp)
    points = [copy(G, g) for g in cp]
    for r in 1:n-1
        for i in 1:n-r
            X = log(G, points[i], points[i+1])
            exp!(G, points[i], points[i], t * X)
        end
    end
    return points[1]
end

R(α) = [cos(α) -sin(α); sin(α) cos(α)]
α(R) = atan(R[2, 1], R[1, 1])
SO2 = SpecialOrthogonalGroup(2)
R2 = TranslationGroup(2)
SE2 = SpecialEuclideanGroup(2)
SO2xR2 = SO2 × R2

cp = [
    ArrayPartition(R(π/2), [0.0, 0.0]),
    ArrayPartition(R(π/4), [0.33, 0.5]),
    ArrayPartition(R(3π/4), [0.66, -0.5]),
    ArrayPartition(R(π), [1.0, 0.0]),
]

# Start a figure and draw the control points
fig = Figure(size = (600, 600))
ax = Axis(fig[1, 1])
positions = [p.x[2] for p in cp]
angles = [α(p.x[1]) for p in cp]
directions = [[cos(a), sin(a)] for a in angles]
scatter!(ax, [p[1] for p in positions], [p[2] for p in positions], color=:blue, markersize=5)
arrows2d!(ax,
    [p[1] for p in positions],
    [p[2] for p in positions],
    [d[1] for d in directions],
    [d[2] for d in directions];
    lengthscale=0.15, color=:blue
)

# Evaluate the Beziér curve at multiple t values
ts1 = range(0, 1, length=35)
ts1_dense = range(0, 1, length=500)
bezier1_points = [deCasteljauLieGroup(SE2, cp, t) for t in ts1]
bezier1_dense_points = [deCasteljauLieGroup(SE2, cp, t) for t in ts1_dense]
bezier1_positions = [p.x[2] for p in bezier1_points]
bezier1_dense_positions = [p.x[2] for p in bezier1_dense_points]
bezier1_angles = [α(p.x[1]) for p in bezier1_points]
bezier1_directions = [[cos(a), sin(a)] for a in bezier1_angles]

lines!(ax, [p[1] for p in bezier1_dense_positions], [p[2] for p in bezier1_dense_positions], color=:red)
scatter!(ax, [p[1] for p in bezier1_positions], [p[2] for p in bezier1_positions], color=:red, markersize=5)
arrows2d!(ax,
    [p[1] for p in bezier1_positions],
    [p[2] for p in bezier1_positions],
    [d[1] for d in bezier1_directions],
    [d[2] for d in bezier1_directions];
    lengthscale=0.15, color=:red, shaftwidth=2.0
)

# Evaluate the Beziér curve at multiple t values
ts1 = range(0, 1, length=35)
ts1_dense = range(0, 1, length=500)
bezier2_points = [deCasteljauLieGroup(SO2xR2, cp, t) for t in ts1]
bezier2_dense_points = [deCasteljauLieGroup(SO2xR2, cp, t) for t in ts1_dense]
bezier2_positions = [p.x[2] for p in bezier2_points]
bezier2_dense_positions = [p.x[2] for p in bezier2_dense_points]
bezier2_angles = [α(p.x[1]) for p in bezier2_points]
bezier2_directions = [[cos(a), sin(a)] for a in bezier2_angles]

lines!(ax, [p[1] for p in bezier2_dense_positions], [p[2] for p in bezier2_dense_positions], color=:green)
scatter!(ax, [p[1] for p in bezier2_positions], [p[2] for p in bezier2_positions], color=:green, markersize=5)
arrows2d!(ax,
    [p[1] for p in bezier2_positions],
    [p[2] for p in bezier2_positions],
    [d[1] for d in bezier2_directions],
    [d[2] for d in bezier2_directions];
    lengthscale=0.15, color=:green, shaftwidth=2.0
)
fig

# export CSV, x,y, u ,v and try to quiver these in TikZ