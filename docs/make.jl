using Documenter
using SkeelBerzins, DifferentialEquations

makedocs(
    modules = [SkeelBerzins, Base.get_extension(SkeelBerzins, :SkeelBerzinsDiffEq)],
    sitename = "SkeelBerzins.jl",
    doctest = false, clean = true,
    warnonly = true,
    format = Documenter.HTML(),
    authors = "Grégoire Pourtier",
    pages = ["Home" => "index.md",
             "Manual" => ["Problem Definition" => "problem_definition.md",
                          "Solvers" => "solvers.md",
                          "Achieve performance" => "performance.md",
                          "Public and Private APIs" => "public_private_APIs.md"
                         ],
             "Examples" => ["101: Linear Diffusion equation" =>"examples/example101.md",
                            "102: Nonlinear Diffusion equation" =>"examples/example102.md",
                            "103: Linear Diffusion equation in cylindrical coordinates" =>"examples/example103.md",
                            "104: Poisson equation" =>"examples/example104.md",
                            "105: Stationary nonlinear diffusion equation" =>"examples/example105.md",
                            "106: System of Reaction-Diffusion equations" =>"examples/example106.md",
                            "107: Linear Diffusion equation in spherical coordinates" =>"examples/example107.md",
                            "201: Interpolation of Partial Derivatives" =>"examples/example201.md",
                            "301: PDE Constrained Optimization" =>"examples/example301.md"
                            ]
            ]
)

deploydocs(
    repo = "github.com/gregoirepourtier/SkeelBerzins.jl.git"
)
