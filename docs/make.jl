#!/usr/bin/env julia
#
#

if "--help" ∈ ARGS
    println(
        """
docs/make.jl

Render the `LieGroups.jl` documentation with optional arguments

Arguments
* `--exclude-tutorials` - exclude the tutorials from the menu of Documenter,
  this can be used if you do not have Quarto installed to still be able to render the docs
  locally on this machine. This option should not be set on CI.
* `--help`              - print this help and exit without rendering the documentation
* `--prettyurls`        – toggle the prettyurls part to true (which is otherwise only true on CI)
* `--quarto`            – run the Quarto notebooks from the `tutorials/` folder before generating the documentation
  this has to be run locally at least once for the `tutorials/*.md` files to exist that are included in
  the documentation (see `--exclude-tutorials`) for the alternative.
  If they are generated once they are cached accordingly.
  Then you can spare time in the rendering by not passing this argument.
  If quarto is not run, some tutorials are generated as empty files, since they
  are referenced from within the documentation. This is currently `getstarted.md`.
""",
    )
    exit(0)
end

run_quarto = "--quarto" in ARGS
run_on_CI = (get(ENV, "CI", nothing) == "true")
tutorials_in_menu = !("--exclude-tutorials" ∈ ARGS)
## Build tutorials menu
tutorials_menu =
    "How to..." => [
        "🚀 Get Started with LieGroups.jl" => "tutorials/getstarted.md",
        "Transition from `GroupManifolds`" => "tutorials/transition.md",
    ]
all_tutorials_exist = true
for (name, file) in tutorials_menu.second
    fn = joinpath(@__DIR__, "src/", file)
    if !isfile(fn) || filesize(fn) == 0 # nonexistent or empty file
        global all_tutorials_exist = false
        if !run_quarto
            @warn "Tutorial $name does not exist at $fn."
            if (!isfile(fn)) && (endswith(file, "getstarted.md"))
                @warn "Generating empty file, since this tutorial is linked to from the documentation."
                touch(fn)
            end
        end
    end
end
if !all_tutorials_exist && !run_quarto && !run_on_CI
    @warn """
        Not all tutorials exist. Run `make.jl --quarto` to generate them. For this run they are excluded from the menu.
    """
    tutorials_in_menu = false
end
if !tutorials_in_menu
    @warn """
    You are either explicitly or implicitly excluding the tutorials from the documentation.
    You will not be able to see their menu entries nor their rendered pages.
    """
    run_on_CI &&
        (@error "On CI, the tutorials have to be either rendered with Quarto or be cached.")
end
#
# (b) if docs is not the current active environment, switch to it
# (from https://github.com/JuliaIO/HDF5.jl/pull/1020/) 
if Base.active_project() != joinpath(@__DIR__, "Project.toml")
    using Pkg
    Pkg.activate(@__DIR__)
    # local temp hack - load ManifoldsBase and Manifold in dev as well
    Pkg.develop(PackageSpec(; path=(@__DIR__) * "/../"))
    Pkg.resolve()
    Pkg.instantiate()
end

# (c) If quarto is set, or we are on CI, run quarto
if run_quarto || run_on_CI
    using CondaPkg
    CondaPkg.withenv() do
        @info "Rendering Quarto"
        tutorials_folder = (@__DIR__) * "/../tutorials"
        # instantiate the tutorials environment if necessary
        Pkg.activate(tutorials_folder)
        # For a breaking release -> also set the tutorials folder to the most recent version
        Pkg.develop(PackageSpec(; path=(@__DIR__) * "/../"))
        Pkg.resolve()
        Pkg.instantiate()
        Pkg.build("IJulia") # build `IJulia` to the right version.
        Pkg.activate(@__DIR__) # but return to the docs one before
        run(`quarto render $(tutorials_folder)`)
        return nothing
    end
else
    # fallback to at least create empty files for tutorials that are directly linked from the docs
    touch(joinpath(@__DIR__, "src/tutorials/getstarted.md"))
end

# (d) load necessary packages for the docs
using Documenter
using DocumenterCitations, DocumenterInterLinks
using LieGroups

# (e) add contributing.md to docs
generated_path = joinpath(@__DIR__, "src")
base_url = "https://github.com/JuliaManifolds/LieGroups.jl/blob/main/"
isdir(generated_path) || mkdir(generated_path)
for (md_file, doc_file) in [("CONTRIBUTING.md", "contributing.md"), ("NEWS.md", "news.md")]
    open(joinpath(generated_path, doc_file), "w") do io
        # Point to source license file
        println(
            io,
            """
            ```@meta
            EditURL = "$(base_url)$(md_file)"
            ```
            """,
        )
        # Write the contents out below the meta block
        for line in eachline(joinpath(dirname(@__DIR__), md_file))
            println(io, line)
        end
    end
end

# (f) finally make docs
bib = CitationBibliography(joinpath(@__DIR__, "src", "references.bib"); style=:alpha)
links = InterLinks(
    "ManifoldsBase" => ("https://juliamanifolds.github.io/ManifoldsBase.jl/stable/"),
    "Manifolds" => ("https://juliamanifolds.github.io/Manifolds.jl/stable/"),
)
makedocs(;
    format=Documenter.HTML(;
        prettyurls=(get(ENV, "CI", nothing) == "true") || ("--prettyurls" ∈ ARGS),
        assets=["assets/favicon.ico", "assets/citations.css", "assets/link-icons.css"],
    ),
    modules=[LieGroups],
    authors="Seth Axen, Mateusz Baran, Ronny Bergmann, Olivier Verdier, and contributors",
    sitename="LieGroups.jl",
    pages=[
        "Home" => "index.md",
        "About" => "about.md",
        (tutorials_in_menu ? [tutorials_menu] : [])...,
        "Lie groups" => [
            "List of Lie groups" => "groups/index.md",
            "General Linear" => "groups/general_linear.md",
            "Heisenberg" => "groups/heisenberg_group.md",
            "Orthogonal group" => "groups/orthogonal_group.md",
            "Power group" => "groups/power_group.md",
            "Product group" => "groups/product_group.md",
            "Semidirect product group" => "groups/semidirect_product_group.md",
            "Special Euclidean group" => "groups/special_euclidean_group.md",
            "Special orthogonal group" => "groups/special_orthogonal_group.md",
            "Special unitary group" => "groups/special_unitary_group.md",
            "Translation group" => "groups/translation_group.md",
            "Unitary group" => "groups/unitary_group.md",
        ],
        "Interfaces" => [
            "Lie group" => "interface/group.md",
            "Lie algebra" => "interface/algebra.md",
            "Group operation" => "interface/operations.md",
            "Group action" => "interface/actions.md",
        ],
        "Contributing to LieGroups.jl" => "contributing.md",
        "Notation" => "notation.md",
        "Changelog" => "news.md",
        "References" => "references.md",
    ],
    plugins=[bib, links],
)
deploydocs(; repo="github.com/JuliaManifolds/LieGroups.jl", push_preview=true)
#back to main env
Pkg.activate()
