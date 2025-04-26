using DLProteinFormats
using Documenter

DocMeta.setdocmeta!(DLProteinFormats, :DocTestSetup, :(using DLProteinFormats); recursive=true)

makedocs(;
    modules=[DLProteinFormats],
    authors="murrellb <murrellb@gmail.com> and contributors",
    sitename="DLProteinFormats.jl",
    format=Documenter.HTML(;
        canonical="https://MurrellGroup.github.io/DLProteinFormats.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/MurrellGroup/DLProteinFormats.jl",
    devbranch="main",
)
