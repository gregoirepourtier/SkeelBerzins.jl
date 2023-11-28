using Documenter
using SkeelBerzins, DifferentialEquations
using ExampleJuggler


function create_doc()
    ExampleJuggler.verbose!(true)

    cleanexamples()

    exampledir = joinpath(@__DIR__, "..", "examples")

    modules = filter(ex -> splitext(ex)[2] == ".jl", basename.(readdir(exampledir)))
    module_examples = @docmodules(exampledir, modules)

    module_examples_reformated = Array{Pair{String,String}}(undef, length(module_examples))
    for i in eachindex(module_examples)
        mod = module_examples[i]
        tmp_new = tmp = replace(replace(first(mod), "_" => ": " ; count=1), "Example" => "")
        cpt_spaces = 0
        for i in 2:length(tmp[2:end-1])
            if isuppercase(tmp[i]) && islowercase(tmp[i+1])
                tmp_new = tmp_new[1:i-1+cpt_spaces]*" "*tmp[i:end]
                cpt_spaces += 1
            end
        end
        module_examples_reformated[i] = Pair(tmp_new, mod[2]) 
    end

    makedocs(; modules = [SkeelBerzins, Base.get_extension(SkeelBerzins, :SkeelBerzinsDiffEq)],
             sitename = "SkeelBerzins.jl",
             checkdocs = :all,
             doctest = false, clean = false,
             warnonly = true,
             format = Documenter.HTML(; mathengine = MathJax3(), 
                                repolink = "https://github.com/gregoirepourtier/SkeelBerzins.jl"),
             authors = "G. Pourtier",
             repo = "https://github.com/gregoirepourtier/SkeelBerzins.jl",
             pages = ["Home" => "index.md",
                      "Manual" => ["Problem Definition" => "problem_definition.md",
                                   "Solvers" => "solvers.md",
                                   "Achieve performance" => "performance.md",
                                   "Public and Private APIs" => "public_private_APIs.md"
                                  ],
                      "Examples" => module_examples_reformated
                    ])

    cleanexamples()

    deploydocs(
        repo = "github.com/gregoirepourtier/SkeelBerzins.jl.git"
    )

end

create_doc()
