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
using Pkg
if Base.active_project() != joinpath(@__DIR__, "Project.toml")
    Pkg.activate(@__DIR__)
end
# load LieGroups in dev mode
Pkg.develop(PackageSpec(; path = (@__DIR__) * "/../"))
Pkg.resolve()
Pkg.instantiate()

# (c) If quarto is set, or we are on CI, run quarto
if run_quarto || run_on_CI
    @info "Rendering Quarto"
    tutorials_folder = (@__DIR__) * "/../tutorials"
    # instantiate the tutorials environment if necessary
    Pkg.activate(tutorials_folder)
    # For a breaking release -> also set the tutorials folder to the most recent version
    Pkg.develop(PackageSpec(; path = (@__DIR__) * "/../"))
    Pkg.resolve()
    Pkg.instantiate()
    Pkg.activate(@__DIR__) # but return to the docs one before
    run(`quarto render $(tutorials_folder)`)
else
    # fallback to at least create empty files for tutorials that are directly linked from the docs
    touch(joinpath(@__DIR__, "src/tutorials/getstarted.md"))
end

# (d) load necessary packages for the docs
using Documenter
using DocumenterCitations, DocumenterInterLinks
using LinearAlgebra
using LieGroups

# (e) add contributing.md and changelog.md to the docs – and link to releases and issues

function add_links(
        line::String, url::String = "https://github.com/JuliaManifolds/LieGroups.jl"
    )
    # replace issues (#XXXX) -> ([#XXXX](url/issue/XXXX))
    while (m = match(r"\(\#([0-9]+)\)", line)) !== nothing
        id = m.captures[1]
        line = replace(line, m.match => "([#$id]($url/issues/$id))")
    end
    # replace ## [X.Y.Z] -> with a link to the release [X.Y.Z](url/releases/tag/vX.Y.Z)
    while (m = match(r"\#\# \[([0-9]+.[0-9]+.[0-9]+)\] (.*)", line)) !== nothing
        tag = m.captures[1]
        date = m.captures[2]
        line = replace(line, m.match => "## [$tag]($url/releases/tag/v$tag) ($date)")
    end
    return line
end

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
            println(io, add_links(line))
        end
    end
end

# (f) finally make docs
bib = CitationBibliography(joinpath(@__DIR__, "src", "references.bib"); style = :alpha)
links = InterLinks(
    "ManifoldsBase" => ("https://juliamanifolds.github.io/ManifoldsBase.jl/stable/"),
    "Manifolds" => ("https://juliamanifolds.github.io/Manifolds.jl/stable/"),
)
makedocs(;
    format = Documenter.HTML(;
        prettyurls = (get(ENV, "CI", nothing) == "true") || ("--prettyurls" ∈ ARGS),
        assets = ["assets/favicon.ico", "assets/citations.css", "assets/link-icons.css"],
        size_threshold_warn = 200 * 2^10, # raise slightly from 100 to 200 KiB
        size_threshold = 300 * 2^10,      # raise slightly 200 to to 300 KiB
    ),
    modules = [LieGroups],
    authors = "Seth Axen, Mateusz Baran, Ronny Bergmann, Olivier Verdier, and contributors",
    sitename = "LieGroups.jl",
    pages = [
        "Home" => "index.md",
        "About" => "about.md",
        (tutorials_in_menu ? [tutorials_menu] : [])...,
        "Lie groups" => [
            "List of Lie groups" => "groups/index.md",
            "Circle Group" => "groups/circle_group.md",
            "General linear group" => "groups/general_linear.md",
            "Heisenberg group" => "groups/heisenberg_group.md",
            "Orthogonal group" => "groups/orthogonal_group.md",
            "Power group" => "groups/power_group.md",
            "Product group" => "groups/product_group.md",
            "Semidirect product group" => "groups/semidirect_product_group.md",
            "Special Euclidean group" => "groups/special_euclidean_group.md",
            "Special linear" => "groups/special_linear.md",
            "Special orthogonal group" => "groups/special_orthogonal_group.md",
            "Special unitary group" => "groups/special_unitary_group.md",
            "Symplectic group" => "groups/symplectic_group.md",
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
    plugins = [bib, links],
)
deploydocs(; repo = "github.com/JuliaManifolds/LieGroups.jl", push_preview = true)
#back to main env
Pkg.activate()
