# Independent series ground truth for `jacobian_exp`, built purely from LieGroups
# primitives (`lie_bracket`, `hat`, `vee`) so the tests validate the closed forms
# against math rather than hard-coded magic numbers. Shared by the SE and SGal tests.

using LieGroups, ManifoldsBase
using LinearAlgebra

if !isdefined(@__MODULE__, :_jacobian_exp_series)

    # small adjoint matrix on the lie algebra X (ad_X), assembled column by column from the lie_bracket,
    # ad_X e_j = [X, E_j]. The input is re-`hat`ed from its coordinates so the bracket sees a
    # single, representation-neutral tangent type.
    function _adjoint_algebra_matrix(G, X)
        𝔤 = LieAlgebra(G)
        n = manifold_dimension(G)
        Xc = hat(𝔤, vee(𝔤, X))
        A = zeros(n, n)
        for j in 1:n
            e_j = zeros(n)
            e_j[j] = 1.0
            A[:, j] = vee(𝔤, lie_bracket(𝔤, Xc, hat(𝔤, e_j)))
        end
        return A
    end

    # Left-trivialized Jacobian of exp via its series J_exp(X) = ∑_k (-ad_X)^k / (k+1)!
    # [Kelly:2025, eq. (28)], an independent check for the closed-form `jacobian_exp`
    function _jacobian_exp_series(G, X; order = 20)
        ad = _adjoint_algebra_matrix(G, X)
        n = size(ad, 1)
        J = Matrix{Float64}(LinearAlgebra.I, n, n)
        term = copy(J)
        for k in 1:order
            term = term * (-ad) ./ (k + 1)  # tₖ = (-ad)ᵏ / (k+1)!
            J .+= term
        end
        return J
    end

end
