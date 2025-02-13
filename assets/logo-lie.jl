using Manifolds, LinearAlgebra, PGFPlotsX, Colors, Distributions, Contour, Random

#
# Settings
#
dark_mode = false

line_offset_brightness = 0.25
patch_opacity = 1.0
geo_opacity = dark_mode ? 1.0 : 1.0
geo_line_width = 20
mesh_line_width = 5
mesh_opacity = dark_mode ? 0.5 : 0.7
logo_colors = [(77, 100, 174), (57, 151, 79), (202, 60, 50), (146, 89, 163)] # Julia colors

rgb_logo_colors = map(x -> RGB(x ./ 255...), logo_colors)
rgb_logo_colors_bright = map(
    x -> RGB((1 + line_offset_brightness) .* x ./ 255...), logo_colors
)
rgb_logo_colors_dark = map(
    x -> RGB((1 - line_offset_brightness) .* x ./ 255...), logo_colors
)

out_file_prefix = dark_mode ? "logo-dark" : "logo"
out_file_ext = ".svg"

#
# Helping functions
#
polar_to_cart(r, θ) = (r * cos(θ), r * sin(θ))

cart_to_polar(x, y) = (hypot(x, y), atan(y, x))

normal_coord_to_vector(M, x, rθ, B) = get_vector(M, x, collect(polar_to_cart(rθ...)), B)

normal_coord_to_point(M, x, rθ, B) = exp(M, x, normal_coord_to_vector(M, x, rθ, B))

function plot_normal_coord!(ax, M, x, B, rs, θs; ncirc=9, options=Dict(), kwargs...)
    for r in rs[2:(end - 1)]
        push!(
            ax,
            Plot3(
                options,
                Coordinates(map(θ -> Tuple(normal_coord_to_point(M, x, [r, θ], B)), θs)),
            ),
        )
    end
    for θ in range(0, 2π; length=ncirc)
        push!(
            ax,
            Plot3(
                options,
                Coordinates(map(r -> Tuple(normal_coord_to_point(M, x, [r, θ], B)), rs)),
            ),
        )
    end
    return ax
end

function plot_patch!(ax, M, x, B, r, θs; options=Dict())
    push!(
        ax,
        Plot3(
            options,
            Coordinates(map(θ -> Tuple(normal_coord_to_point(M, x, [r, θ], B)), θs)),
        ),
    )
    return ax
end

# skip before end removed end-(m-1)...end-1, so m elements
function plot_geodesic!(ax, M, x, y; n=100, m=0, options=Dict())
    γ = shortest_geodesic(M, x, y)
    T = range(0, 1; length=n)
    pts = γ.(T)
    cut = [pts[1:(end - (m + 1))]..., last(pts)]
    push!(ax, Plot3(options, Coordinates(Tuple.(cut))))
    return ax
end

#
# Prepare document
#
resize!(PGFPlotsX.CUSTOM_PREAMBLE, 0)
push!(PGFPlotsX.CUSTOM_PREAMBLE, raw"\pgfplotsset{scale=6.0}")
push!(PGFPlotsX.CUSTOM_PREAMBLE, raw"\usetikzlibrary{bending}")
push!(PGFPlotsX.CUSTOM_PREAMBLE, raw"\usetikzlibrary{arrows.meta}")
push!(PGFPlotsX.CUSTOM_PREAMBLE, raw"\pgfplotsset{roundcaps/.style={line cap=round}}")
push!(
    PGFPlotsX.CUSTOM_PREAMBLE,
    raw"\pgfplotsset{meshlinestyle/.style={dash pattern=on 0pt off 2.5\pgflinewidth, line cap=round}}",
)
if dark_mode
    push!(PGFPlotsX.CUSTOM_PREAMBLE, raw"\pagecolor{black}")
end
S = Sphere(2)

center = normalize([1, 1, 1])
x, y, z = eachrow(Matrix{Float64}(I, 3, 3))
γ1 = shortest_geodesic(S, center, z)
γ2 = shortest_geodesic(S, center, x)
γ3 = shortest_geodesic(S, center, y)
p1 = γ1(1)
p2 = γ2(1)
p3 = γ3(1)

#
# Setup Axes
if dark_mode
    tp = @pgf Axis({
        axis_lines = "none",
        axis_equal,
        view = "{135}{35}",
        zmin = -0.05,
        zmax = 1.0,
        xmin = 0.0,
        xmax = 1.0,
        ymin = 0.0,
        ymax = 1.0,
    })
else
    tp = @pgf Axis({
        axis_lines = "none",
        axis_equal,
        view = "{135}{35}",
        zmin = -0.05,
        zmax = 1.0,
        xmin = 0.0,
        xmax = 1.0,
        ymin = 0.0,
        ymax = 1.0,
    })
end
rs = range(0, π / 5; length=6)
θs = range(0, 2π; length=100)

#
# Plot manifold patches
# rotate color indices
indices = [3, 4, 2] # original [2,3,4]
patch_colors = rgb_logo_colors[indices]
patch_colors_line =
    dark_mode ? rgb_logo_colors_bright[indices] : rgb_logo_colors_dark[indices]
base_points = [p1, p2, p3]
basis_vectors = [log(S, p1, p2), log(S, p2, p1), log(S, p3, p1)]
for i in eachindex(base_points)
    xi = base_points[i]
    B = DiagonalizingOrthonormalBasis(basis_vectors[i])
    basis = get_basis(S, xi, B)
    optionsP = @pgf {fill = patch_colors[i], draw = "none", opacity = patch_opacity}
    plot_patch!(tp, S, xi, basis, π / 5, θs; options=optionsP)
    optionsL = @pgf {
        meshlinestyle,
        color = dark_mode ? "white" : "black",
        line_width = mesh_line_width,
        opacity = mesh_opacity,
    }
    plot_normal_coord!(tp, S, xi, basis, rs, θs; options=optionsL)
end

#
# Plot geodesics
options = @pgf {
    opacity = geo_opacity,
    no_markers,
    roundcaps,
    line_width = geo_line_width,
    color = dark_mode ? "white" : "black",
    "-{Stealth[length=6.195cm,sharp,bend=-9]}", #length=8mm,width=2mm,inset=3mm
}
plot_geodesic!(tp, S, base_points[1], base_points[2]; m=25, options=options)
plot_geodesic!(tp, S, base_points[3], base_points[1]; m=25, options=options)
plot_geodesic!(tp, S, base_points[2], base_points[3]; m=25, options=options)

#=
push!(
    tp,
    raw"\node [scale=50, color=" *
    "$(dark_mode ? "white" : "black")] at (0,0,0) " *
    raw"{$\circ$};",
)
=#

#
# Export Logo.
out_file = "$(out_file_prefix)$(out_file_ext)"
pgfsave(out_file, tp)
pgfsave("$(out_file_prefix).pdf", tp)
