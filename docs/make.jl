using LieGroups
using Documenter

DocMeta.setdocmeta!(LieGroups, :DocTestSetup, :(using LieGroups); recursive=true)

makedocs(;
    modules=[LieGroups],
    authors="Yueh-Hua Tu",
    repo="https://github.com/yuehhua/LieGroups.jl/blob/{commit}{path}#{line}",
    sitename="LieGroups.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://yuehhua.github.io/LieGroups.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/yuehhua/LieGroups.jl",
    devbranch="main",
)
