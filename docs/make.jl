using LieGroup
using Documenter

DocMeta.setdocmeta!(LieGroup, :DocTestSetup, :(using LieGroup); recursive=true)

makedocs(;
    modules=[LieGroup],
    authors="Yueh-Hua Tu",
    repo="https://github.com/yuehhua/LieGroup.jl/blob/{commit}{path}#{line}",
    sitename="LieGroup.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://yuehhua.github.io/LieGroup.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/yuehhua/LieGroup.jl",
    devbranch="main",
)
