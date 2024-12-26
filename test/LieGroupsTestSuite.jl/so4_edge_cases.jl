𝔰𝔬4_edges_cases_explog = [
    [0, 0, π, 0, 0, π],  # θ = (π, π)
    [0, 0, π, 0, 0, 0],  # θ = (π, 0)
#    [0, 0, π / 2, 0, 0, π],  # θ = (π, π/2) TODO: Check why this is different from the previous version in Manifolds.jl
    [0, 0, π / 2, 0, 0, 0],  # θ = (0, π/2)
    [0, 0, π, 0, 0, 0] ./ 2,  # θ = (π/2, 0)
    [0, 0, π, 0, 0, π] ./ 2,  # θ = (π/2, π/2)
    [0, 0, 0, 0, 0, 0],  # θ = (0, 0)
    [0, 0, 1, 0, 0, 1] .* 1e-100, # α = β ≈ 0
    [0, 0, 1, 0, 0, 1] .* 1e-6, # α = β ⩰ 0
    [0, 0, 10, 0, 0, 1] .* 1e-6, # α ⪆ β ⩰ 0
    [0, 0, π / 4, 0, 0, π / 4 - 1e-6], # α ⪆ β > 0
]
