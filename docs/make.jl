using Documenter
using SkeelBerzins

makedocs(
    modules = [SkeelBerzins],
    sitename = "SkeelBerzins.jl",
    doctest = false, clean = true,
    format = Documenter.HTML(),
    authors = "Grégoire Pourtier",
    pages = ["Home" => "index.md",
             "Manual" => ["Problem Definition" => "problem_definition.md",
                          "Solvers" => "solvers.md",
                          "Achieve performance" => "performance.md",
                          "Private API" => "private_API.md"
                         ],
             "Examples" => "examples.md"
            ]
)

deploydocs(
    repo = "github.com/gregoirepourtier/SkeelBerzins.jl.git"
)
