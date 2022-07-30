# Reference: https://github.com/artivis/manif/blob/devel/examples/se3_localization.py

using LieGroups
using LinearAlgebra
using Plots

const NUMBER_OF_LMKS_TO_MEASURE = 5
const STEPS = 10
const LANDMARKS = [2.0  0.0  0.0;
                   3.0 -1.0 -1.0;
                   2.0 -1.0  1.0;
                   2.0  1.0  1.0;
                   2.0  1.0 -1.0]

ϵ(x) = rand(x)

function main()
    X_simulation = identity(SE{3})
    X = identity(SE{3})
    X_unfiltered = identity(SE{3})
    P = zeros(dof(SE{3}), dof(SE{3}))

    u_nom = [0.1, 0.0, 0.0, 0.0, 0.0, 0.05]
    u_σs = 0.1
    U = diagm(repeat([abs2.(u_σs)], dof(SE{3})))

    # Define the beacon's measurements
    measurements = zeros(NUMBER_OF_LMKS_TO_MEASURE, 3)

    y_σs = 0.01
    R = diagm(repeat([abs2.(y_σs)], dim(SE{3})))

    # CONFIGURATION DONE

    @show log(X_simulation)

    # START TEMPORAL LOOP

    for _ in 1:STEPS
        # I. Simulation
        u_noisy = u_nom + u_σs .* ϵ(dof(SE{3}))              # simulate noise

        u_simu = se{3}(u_nom)
        u_est = se{3}(u_noisy)
        u_unfilt = se{3}(u_noisy)

        # first we move
        X_simulation = X_simulation ⊕ u_simu                # overloaded X.rplus(u) = X * exp(u)

        # then we measure all landmarks
        for i in 1:NUMBER_OF_LMKS_TO_MEASURE
            y = inv(X_simulation) ⋉ LANDMARKS[i, :]          # landmark measurement, before adding noise

            # simulate noise
            measurements[i, :] .= y + y_σs .* ϵ(dim(SE{3}))  # measure landmarks with noise
        end

        # II. Estimation

        # First we move
        J_x, J_u = jacobian(⊕, X, u_est)
        X = X ⊕ u_est                                       # X * exp(u), with Jacobians

        P = J_x * P * J_x' + J_u * U * J_u'

        # Then we correct using the measurements of each lmk
        for i in 1:NUMBER_OF_LMKS_TO_MEASURE
            b = LANDMARKS[i, :]                              # lmk coordinates in world frame
            y = measurements[i, :]                           # lmk measurement, noisy

            # expectation
            e = inv(X) ⋉ b                                   # note: e = R.tr * ( b - t ), for X = (R,t).
            J_xi_x = jacobian(inv, X)
            J_e_xi = jacobian(⋉, inv(X), b)
            H = J_e_xi * J_xi_x                              # Jacobian of the measurements wrt the robot pose. note: H = J_e_x = J_e_xi * J_xi_x
            E = H * P * H'

            # innovation
            z = y - e
            Z = E + R

            # Kalman gain
            K = P * H' * inv(Z)                              # K = P * H.tr * ( H * P * H.tr + R).inv

            # Correction step
            dx = K * z                                       # dx is in the tangent space at X

            # Update
            X = X ⊕ se{3}(dx)                               # overloaded X.rplus(dx) = X * exp(dx)
            P = P - K * Z * K'
        end

        # III. Unfiltered

        # move also an unfiltered version for comparison purposes
        X_unfiltered = X_unfiltered ⊕ u_unfilt

        # IV. Results

        @show log(X_simulation)
        @show log(X)
        @show log(X_unfiltered)
        println("========================================")
    end
end

main()

scatter(LANDMARKS[:,1], LANDMARKS[:,2], LANDMARKS[:,3], markercolor=:red)
