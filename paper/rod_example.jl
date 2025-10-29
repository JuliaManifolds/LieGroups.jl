using LinearAlgebra, Plots

using LieGroups, Manifolds, RecursiveArrayTools

# ---------- Physical parameters ----------

ϵ0 = 8.854e-12       # vacuum permittivity (F/m)
L   = 1.0             # rod length (m)
M   = 1.0             # rod mass (kg)
λ   = 1e-6            # linear charge density (C/m)
Icm = M * L^2 / 12    # moment of inertia about center
Nq  = 20              # number of sample points along the rod for integration

# Point charges (Coulombs)

q1, q2 = +1e-6, -1e-6
R1 = [-1.0, 0.0]      # position of +q1
R2 = [ 1.0, 0.0]      # position of -q2

# ---------- Electric potential and field ----------

function E_field(pos)
    r1 = pos .- R1
    r2 = pos .- R2
    E = (q1 .* r1 ./ norm(r1)^3 .+ q2 .* r2 ./ norm(r2)^3) ./ (4π*ϵ0)
    return E
end

# ---------- Force and torque on the rod ----------

function force_and_torque(r, Fθ)
    u, u⊥ = eachrow(Fθ)
    svals = range(-L/2, L/2; length=Nq)
    F = zeros(2)
    τ = 0.0

    for s in svals
        pos = r .+ s .* u
        E = E_field(pos)
        F .+= λ * E * (L/(Nq-1))           # integrate via trapezoidal rule spacing
        τ += λ * s * dot(u⊥, E) * (L/(Nq-1))
    end
    return F, τ
end

# ---------- Equations of motion ----------

function dynamics(state, t)
    SE2 = SpecialEuclideanGroup(2)
    se2 = LieAlgebra(SE2)
    p, dp = state.x
    Fθ, r = p.x
    F, τ = force_and_torque(r, Fθ)
    ax, ay = F ./ M
    α = τ / Icm
    ddx = hat(se2, [α, ax, ay], ArrayPartition)
    return ArrayPartition(dp, ddx)
end

# ---------- Time integration (Heun) ----------

function integrator_step!(A::GroupAction, y, f, t, dt)
    Ie = Identity(A.group)
    F1 = (dt / 2) * f(y, t)

    tmp = zero_vector(A.manifold, y)
    diff_group_apply!(A, tmp, Ie, y, F1)
    y2 = exp(A.manifold, y, tmp)
    F2 = dt * f(y2, t + dt/2)
    diff_group_apply!(A, tmp, Ie, y, F2)
    # @show y.x[1].x[1]
    # @show tmp.x[1].x[1]
    exp!(base_manifold(A), y, y, tmp)
    # @show y.x[1].x[1]
    println()
end

function simulate_trajectory(angle, x0, rot_vel, v0)
    dt = 0.001
    tmax = 10.0
    steps = Int(round(tmax/dt))

    SE2 = SpecialEuclideanGroup(2)
    se2 = LieAlgebra(SE2)
    G = ProductLieGroup(SE2, LieGroup(se2, AdditionGroupOperation()))
    action = GroupAction(G, G, LeftGroupOperationAction())

    ISE2 = identity_element(SE2, ArrayPartition)
    # initial position and angle
    state_SE2 = exp(SE2.manifold, ISE2, hat(se2, [angle, x0...], ArrayPartition))
    initial_velocity = hat(se2, [rot_vel, v0...], ArrayPartition)
    state = ArrayPartition(state_SE2, initial_velocity)

    traj = typeof(state)[]
    for i in 1:steps
        push!(traj, copy(G, state))
        integrator_step!(action, state, dynamics, (i-1)*dt, dt)
    end

    centers = [traj[i].x[1].x[2] for i in eachindex(traj)]
    Fθs = [traj[i].x[1].x[1] for i in eachindex(traj)]
    return centers, Fθs
end

# ---------- Simulation ----------
function sim()
    angles = [0.1, 0.3, 0.6]
    positions = [[0.0, -1.5], [0.0, -1.0], [0.5, -1.5]]

    plt = plot(title="Charged rod trajectories near ±q charges",
               xlabel="x (m)", ylabel="y (m)", aspect_ratio=1, legend=false)
    
    for (ang, pos) in zip(angles, positions)
        centers, Fθs = simulate_trajectory(ang, pos, 0.0, [0.0, 0.0])
        x = [c[1] for c in centers]
        y = [c[2] for c in centers]
        plot!(plt, x, y, label="θ₀=$(ang), pos=$(pos)")

        # Draw rods at 10 roughly equidistant times
        nsteps = length(centers)
        idxs = round.(Int, range(1, nsteps, length=10))
        for i in idxs
            c = centers[i]
            u = Fθs[i][1, :]           # rod's local x-axis
            p1 = c .- 0.5L .* u
            p2 = c .+ 0.5L .* u
            @show u
            plot!(plt, [p1[1], p2[1]], [p1[2], p2[2]], color=:black, lw=1, label=false)
        end
    end

    scatter!(plt, [R1[1] R2[1]], [R1[2] R2[2]],
             marker_z=[q1 q2], color=[:red :blue],
             label=["+q" "-q"], markersize=6)
    display(plt)
end

sim()
