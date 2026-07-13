using LieGroups, ManifoldsBase, Random, Test, RecursiveArrayTools
using LieGroups: SpecialGalileanGroup

using StaticArrays
using LinearAlgebra

# --- independent ground truths built from LieGroups primitives (no magic numbers) ---

# small adjoint matrix ad_X and the right-Jacobian series ground truth for `jacobian_exp`
include("jacobian_exp_series_reference.jl")

# reference Lie bracket via the matrix commutator of the 5×5 screw representation
# (see the `hat` docstring), an independent check for the closed-form `lie_bracket`
_sgal3_skew(a) = [0.0 -a[3] a[2]; a[3] 0.0 -a[1]; -a[2] a[1] 0.0]
_sgal3_screw(c) = [_sgal3_skew(c[7:9]) c[4:6] c[1:3]; zeros(1, 3) 0.0 c[10]; zeros(1, 5)]
_sgal3_coords(Z) = [Z[1:3, 5]; Z[1:3, 4]; Z[3, 2]; Z[1, 3]; Z[2, 1]; Z[4, 5]]
function _lie_bracket_ref(G, X, Y)
    𝔤 = LieAlgebra(G)
    Xm, Ym = _sgal3_screw(vee(𝔤, X)), _sgal3_screw(vee(𝔤, Y))
    return hat(𝔤, _sgal3_coords(Xm * Ym - Ym * Xm))
end

@testset "Special Galilean" begin
    𝔰 = sqrt(2)
    fcts = [
        compose,
        exp,
        get_vector,
        hat,
        identity_element,
        inner,
        inv,
        # is_flat,
        is_identity,
        jacobian_exp,
        lie_bracket,
        log,
        norm,
        rand,
        show,
        vee,
        apply,
    ]

    @testset "SGal(3)" begin
        G3p = SpecialGalileanGroup(3)
        hL1 = ArrayPartition(
            ArrayPartition(
                1 / 𝔰 * [1.0 1.0 0.0; -1.0 1.0 0.0; 0.0 0.0 𝔰], [1 / 𝔰, 0.0, 0.0]
            ),
            ArrayPartition([0.0, 0.0, 0.0], [0.0]),
        )
        hL2 = ArrayPartition(
            ArrayPartition([0.0 -1.0 0.0; 1.0 0.0 0.0; 0.0 0.0 1.0], [0.0, 1.0, 0.0]),
            ArrayPartition([0.0, 0.1, 0.0], [0.1]),
        )
        hL3 = ArrayPartition(
            ArrayPartition([1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0], [1.0, 1.0, 0.0]),
            ArrayPartition([0.0, 0.0, 1.0], [1.0]),
        )
        hL4 = identity_element(G3p, StaticArray)

        YL1 = ArrayPartition(
            ArrayPartition([0.0 -0.23 0.0; 0.23 0.0 0.0; 0.0 0.0 0.0], [0.0, 1.0, 0.0]),
            ArrayPartition([0.0, 0.0, 0.0], [1.0]),
        )
        YL2 = ArrayPartition(
            ArrayPartition([0.0 0.3 0.0; -0.3 0.0 0.0; 0.0 0.0 0.0], [1.0, 1.0, 0.0]),
            ArrayPartition([1.0, 0.3, 0.0], [0.1]),
        )
        YL3 = ArrayPartition(
            ArrayPartition([9.0 0.1 0.0; -0.1 0.0 0.0; 0.0 0.0 0.0], [1.0, 0.0, 0.0]),
            ArrayPartition([0.0, 0.0, 0.3], [0.1]),
        )
        YL4 = zero(hL4)
        hP = [hL1, hL2, hL3, hL4]
        YP = [YL1, YL2, YL3, YL4]
        for G in [G3p], (pts, vec) in zip([hP], [YP])
            properties = Dict(
                :Name => "Special Galilean Group SGal(3)",
                :Points => pts,
                :Vectors => vec,
                :Rng => Random.MersenneTwister(),
                :Functions => fcts,
            )
            expectations = Dict(
                :atol => 1.0e-14,
                # :repr => "SpecialGalileanGroup(3)",
                # :is_flat => false
                # jacobian_exp (on vec[1]) and lie_bracket (on vec[1], vec[2]) are validated
                # against independent ground truths built from LieGroups primitives (see the
                # helpers at the top of this file) rather than hard-coded magic numbers
                :jacobian_exp => _jacobian_exp_series(G, vec[1]),
                :lie_bracket => _lie_bracket_ref(G, vec[1], vec[2]),
            )
            LieGroups.Test.test_lie_group(G, properties, expectations)
        end
    end

    # jacobian_exp: exercise BOTH the truncated-series branch (φ < 0.15) and the
    # closed-form branch, each validated against the series ground truth
    @testset "jacobian_exp series and closed-form branches" begin
        G = SpecialGalileanGroup(3)
        g = identity_element(G)
        mkX = φ -> ArrayPartition(
            ArrayPartition([0.0 -φ 0.0; φ 0.0 0.0; 0.0 0.0 0.0], [1.0, 0.5, 0.0]),
            ArrayPartition([0.3, 0.0, 0.2], [0.4]),
        )
        for φ in (1.0e-3, 0.23)  # 1e-3 → series branch (φ<0.15), 0.23 → closed form
            X = mkX(φ)
            @test isapprox(jacobian_exp(G, g, X), _jacobian_exp_series(G, X); atol = 1.0e-12)
        end
    end

    @testset "Test SGal(3) SArray" begin
        G = SpecialGalileanGroup(3)
        ε = identity_element(G, StaticArray)
        c = SA[1.0, 2.0, 3.0, 0.1, 0.2, 0.3, 0.01, 0.02, 0.03, 0.5]
        X = hat(LieAlgebra(G), c, typeof(ε))
        @test typeof(X) == typeof(ε)
        c1 = vee(LieAlgebra(G), X)
        @test c1 ≈ c
        @test typeof(c1) == typeof(c)
        g = exp(G, X)
        @test typeof(g) == typeof(ε)
        # pack c manually into a matrix to test against the matrix exponential
        Xmat = [
            0.0  -0.03 0.02 0.1 1.0;
            0.03  0.0 -0.01 0.2 2.0;
            -0.02  0.01 0.0  0.3 3.0;
            0.0   0.0  0.0  0.0 0.5;
            0.0   0.0  0.0  0.0 0.0
        ]
        gmat = exp(Xmat)
        @test g.x[1].x[1] ≈ gmat[1:3, 1:3]
        @test g.x[1].x[2] ≈ gmat[1:3, 4]
        @test g.x[2].x[1] ≈ gmat[1:3, 5]
        @test g.x[2].x[2][1] ≈ gmat[4, 5]
    end

    @testset "Test SGal(3) extras" begin
        G = SpecialGalileanGroup(3)
        ε = identity_element(G, StaticArray)
        c = SA[1.0, 2.0, 3.0, 0.1, 0.2, 0.3, 0.01, 0.02, 0.03, 0.5]
        X1 = hat(LieAlgebra(G), c, typeof(ε))
        g1 = exp(G, X1)
        inv_g1 = inv(G, g1)
        ε1 = compose(G, inv_g1, g1)

        ε = identity_element(G)
        c = [1.0, 2.0, 3.0, 0.1, 0.2, 0.3, 0.01, 0.02, 0.03, 0.5]
        X2 = hat(LieAlgebra(G), c, typeof(ε))
        g2 = identity_element(G)
        exp!(G, g2, X2)
        inv_g2 = identity_element(G)
        inv!(G, inv_g2, g2)
        ε2 = identity_element(G)
        compose!(G, ε2, inv_g2, g2)

        @test isapprox(G, X1, X2)
        @test isapprox(G, g1, g2)
        @test isapprox(G, inv_g1, inv_g2)
        @test isapprox(G, ε1, ε, atol = 1.0e-14)
        @test isapprox(G, ε2, ε, atol = 1.0e-14)

    end
end
