using DifferentiableOSQP
using Documenter

DocMeta.setdocmeta!(DifferentiableOSQP, :DocTestSetup, :(using DifferentiableOSQP); recursive=true)

makedocs(;
    modules=[DifferentiableOSQP],
    authors="Devansh Agrawal",
    repo="https://github.com/dev10110/DifferentiableOSQP.jl/blob/{commit}{path}#{line}",
    sitename="DifferentiableOSQP.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://dev10110.github.io/DifferentiableOSQP.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/dev10110/DifferentiableOSQP.jl",
    devbranch="main",
)
