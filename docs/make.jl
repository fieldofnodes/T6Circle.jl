using T6Circle
using Documenter

DocMeta.setdocmeta!(T6Circle, :DocTestSetup, :(using T6Circle); recursive=true)

makedocs(;
    modules=[T6Circle],
    authors="Jonathan Miller jonathan.miller@fieldofnodes.com",
    sitename="T6Circle.jl",
    format=Documenter.HTML(;
        canonical="https://fieldofnodes.github.io/T6Circle.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/fieldofnodes/T6Circle.jl",
    devbranch="main",
)
