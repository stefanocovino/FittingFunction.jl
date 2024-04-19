using FittingFunction
using Documenter

DocMeta.setdocmeta!(FittingFunction, :DocTestSetup, :(using FittingFunction); recursive=true)

makedocs(;
    modules=[FittingFunction],
    authors="Stefano Covino <stefano.covino@inaf.it> and contributors",
    sitename="FittingFunction.jl",
    format=Documenter.HTML(;
        canonical="https://stefanocovino.github.io/FittingFunction.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/stefanocovino/FittingFunction.jl",
    devbranch="main",
)
