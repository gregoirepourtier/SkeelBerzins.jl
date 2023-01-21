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
             "Examples" => ["101: Linear Diffusion equation (DiffEq)" =>"example101.md",
                            "102: Linear Diffusion equation" =>"example102.md",
                            "103: Nonlinear Diffusion equation (DiffEq)" =>"example103.md",
                            "104: Nonlinear Diffusion equation" =>"example104.md",
                            "105: Linear Diffusion equation in cylindrical coordinates (DiffEq)" =>"example105.md",
                            "106: Linear Diffusion equation cylindrical coordinates" =>"example106.md",
                            "107: Poisson equation" =>"example107.md",
                            "108: Stationary nonlinear diffusion equation" =>"example108.md",
                            "109: System of Reaction-Diffusion equations (DiffEq)" =>"example109.md",
                            "110: System of Reaction-Diffusion equations" =>"example110.md",
                            "111: Linear Diffusion equation in spherical coordinates (DiffEq)" =>"example111.md",
                            "112: Linear Diffusion equation in spherical coordinates" =>"example112.md",
                            ]
            ]
)

deploydocs(
    repo = "github.com/gregoirepourtier/SkeelBerzins.jl.git"
)
