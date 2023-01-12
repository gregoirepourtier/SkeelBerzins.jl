using Documenter
using SkeelBerzins, SciMLBase

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
                          "Public and Private APIs" => "public_private_APIs.md"
                         ],
             "Examples" => "examples.md"
            ]
)

deploydocs(
    repo = "github.com/gregoirepourtier/SkeelBerzins.jl.git"
)
