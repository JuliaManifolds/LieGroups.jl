using Pkg; Pkg.activate(@__DIR__); Pkg.instantiate()
cd(@__DIR__)
using LinearAlgebra
using LieGroups, Manifolds, RecursiveArrayTools
using CairoMakie
using DataFrames, CSV

# ---------- Physical parameters ----------

Base.@kwdef struct SystemParameters
    group = ProductLieGroup(SpecialOrthogonalGroup(2), TranslationGroup(2))
    lie_algebra = LieAlgebra(group)
    L = 0.67             # rod length (m)
    M = 1.0             # rod mass (kg)
    λ = 1e-6            # linear charge density (C/m)
    Icm = M * L^2 / 12  # moment of inertia about center
    qs = [+1e-6, -1e-6] # point charges (Coulombs)
    Rs = [[-1.0, 0.0], [1.0, 0.0]] # positions of point charges
    initial_rod_angle = 0.0
    initial_rod_position = [0.0, -1.0]
    initial_rod_angular_velocity = 0.1
    initial_rod_velocity = [0.0, 0.0]
    tmax::Real = 30.0
end

ϵ0 = 8.854e-12  # vacuum permittivity (F/m)
Nq = 20         # number of sample points along the rod for integration


# ---------- Electric potential and field ----------

function E_field(sp::SystemParameters, pos)
    E = zeros(2)
    for i in eachindex(sp.qs)
        ri = pos .- sp.Rs[i]
        E .+= (sp.qs[i] .* ri ./ norm(ri)^3) ./ (4π*ϵ0)
    end
    return E
end

# ---------- Equations of motion ----------
# first block for the paper

function force_and_torque(sp::SystemParameters, Fθ, r)
    u, u⊥ = eachrow(Fθ)
    F, τ = zeros(2), 0.0
    for s in range(-sp.L/2, sp.L/2; length=Nq)
        E = E_field(sp, r .+ s .* u)
        F .+= sp.λ * E * (sp.L/(Nq-1))
        τ += sp.λ * s * dot(u⊥, E) * (sp.L/(Nq-1))
    end
    return F, τ
end

function dynamics(state, sp::SystemParameters, t)
    p, dp = state.x
    F, τ = force_and_torque(sp, p.x...)
    a = F ./ sp.M
    ddx = hat(sp.lie_algebra, [τ / sp.Icm, a...], ArrayPartition)
    return ArrayPartition(dp, ddx)
end

# ---------- Time integration (Heun) ----------
# second block for the paper

function integrator_step!(A::GroupAction,
    sp::SystemParameters, y, f, t, dt)
    Ie = Identity(A.group)
    F1 = (dt / 2) * f(y, sp, t)
    tmp = zero_vector(A.manifold, y)
    diff_group_apply!(A, tmp, Ie, y, F1)
    y2 = exp(A.manifold, y, tmp)
    F2 = dt * f(y2, sp, t + dt/2)
    diff_group_apply!(A, tmp, Ie, y, F2)
    exp!(base_manifold(A), y, y, tmp)
end


# ---------- Trajectory simulation ----------

function simulate_trajectory(BG, sp::SystemParameters, angle, x0, rot_vel, v0)
    dt = 0.001
    steps = Int(round(sp.tmax/dt))

    ba = LieAlgebra(BG)
    G = ProductLieGroup(BG, LieGroup(ba, AdditionGroupOperation()))
    action = GroupAction(G, G, LeftGroupOperationAction())

    IBG = identity_element(BG, ArrayPartition)
    # initial position and angle
    state_BG = exp(BG.manifold, IBG, hat(ba, [angle, x0...], ArrayPartition))
    initial_velocity = hat(ba, [rot_vel, v0...], ArrayPartition)
    state = ArrayPartition(state_BG, initial_velocity)

    traj = typeof(state)[]
    for i in 1:steps
        push!(traj, copy(G, state))
        integrator_step!(action, sp, state, dynamics, (i-1)*dt, dt)
    end

    centers = [traj[i].x[1].x[2] for i in eachindex(traj)]
    Fθs = [traj[i].x[1].x[1] for i in eachindex(traj)]
    return centers, Fθs
end

# ---------- Simulation ----------
function sim()
    fig = Figure(size = (800, 800))
    ax = Axis(fig[1, 1],
        title = "Charged rod trajectories near ±q charges",
        xlabel = "x (m)", ylabel = "y (m)",
        aspect = DataAspect(),
    )

    SE2 = SpecialEuclideanGroup(2)
    SOxT = ProductLieGroup(SpecialOrthogonalGroup(2), TranslationGroup(2))
    sp = SystemParameters()

    sps = [SystemParameters(; initial_rod_angle = 0.1, initial_rod_position=[0.0, -1.5], tmax=135.0),
           SystemParameters(; initial_rod_angle = 0.2, initial_rod_position=[-0.2, +1.0], tmax=30.0, initial_rod_velocity = [-0.12, -0.15]),
           SystemParameters(; initial_rod_angle = 0.6, initial_rod_position=[-0.5, -0.5], tmax=30.0, initial_rod_velocity = [-0.04, 0.12])]

    export_folder = "data/"
    colors = [:red, :blue, :green]

    for (isp, sp) in enumerate(sps)
        ang = sp.initial_rod_angle
        pos = sp.initial_rod_position
        centers, Fθs = simulate_trajectory(SOxT, sp, ang, pos, sp.initial_rod_angular_velocity, sp.initial_rod_velocity)
        x = [c[1] for c in centers]
        y = [c[2] for c in centers]
        lines!(ax, x, y, color = colors[isp])

        # Draw rods at 10 roughly equidistant times
        nsteps = length(centers)
        idxs = round.(Int, range(1, nsteps, length=20))
        for i in idxs
            c = centers[i]
            u = Fθs[i][1, :]
            p1 = c .- 0.5sp.L .* u
            p2 = c .+ 0.5sp.L .* u
            lines!(ax, [p1[1], p2[1]], [p1[2], p2[2]], color=:black, linewidth=1)
        end

        df_dense = DataFrame(x = vcat(x[begin:500:end], x[end]), y = vcat(y[begin:500:end], y[end]))
        df_rod = DataFrame(
            x = x[idxs],
            y = y[idxs],
            up = [0.5sp.L * Fθs[i][1, 1] for i in idxs],
            vp = [0.5sp.L * Fθs[i][1, 2] for i in idxs],
            um = [-0.5sp.L * Fθs[i][1, 1] for i in idxs],
            vm = [-0.5sp.L * Fθs[i][1, 2] for i in idxs],
        )
        CSV.write(export_folder * "ex2_rod_$(isp)_dense.csv", df_dense)
        CSV.write(export_folder * "ex2_rod_$(isp)_sparse.csv", df_rod)
    end

    # scatter the two charges
    scatter!(ax, [sp.Rs[1][1], sp.Rs[2][1]],
                 [sp.Rs[1][2], sp.Rs[2][2]],
                 color = [:red, :blue], markersize = 12)

    display(fig)
    save("img/rod_trajectories_makie.png", fig)
end

sim()