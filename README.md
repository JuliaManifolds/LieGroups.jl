# LieGroups

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://yuehhua.github.io/LieGroups.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://yuehhua.github.io/LieGroups.jl/dev)
[![Build Status](https://github.com/yuehhua/LieGroups.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/yuehhua/LieGroups.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/yuehhua/LieGroups.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/yuehhua/LieGroups.jl)

LieGroups provides Lie group objects and their operations.

## Usage

```julia
julia> using LieGroups

julia> θ = 45/180*π
0.7853981633974483

julia> alg = so{2}([θ])
so{2, Vector{Float64}}([0.7853981633974483])

julia> g = exp(∧(alg))
SO{2}(A=[0.7071067811865476 -0.7071067811865475; 0.7071067811865475 0.7071067811865477])

julia> g ⋉ [1, 0]
2-element Vector{Float64}:
 0.7071067811865476
 0.7071067811865475

julia> (g * g) ⋉ [1, 0]
2-element Vector{Float64}:
 2.220446049250313e-16
 1.0

julia> alg2 = log(g * g)
so{2, Vector{Float64}}([1.5707963267948966])

julia> alg2 == alg + alg
true
```
