using Documenter
using SkeelBerzins, DifferentialEquations

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
             "Examples" => ["101: Linear Diffusion equation" =>"example101.md",
                            "102: Nonlinear Diffusion equation" =>"example102.md",
                            "103: Linear Diffusion equation in cylindrical coordinates" =>"example103.md",
                            "104: Poisson equation" =>"example104.md",
                            "105: Stationary nonlinear diffusion equation" =>"example105.md",
                            "106: System of Reaction-Diffusion equations" =>"example106.md",
                            "107: Linear Diffusion equation in spherical coordinates" =>"example107.md",
                            "201: Interpolation of Partial Derivatives" =>"example201.md",
                            "301: PDE Constrained Optimization" =>"example301.md"
                            ]
            ]
)

deploydocs(
    repo = "github.com/gregoirepourtier/SkeelBerzins.jl.git"
)
