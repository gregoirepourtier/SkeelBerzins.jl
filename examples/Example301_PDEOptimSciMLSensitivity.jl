#=

# Example 301: PDE Optim SciML Sensitivity

This example is used to show the integration of the package in the SciML ecosystem.
For more details on the problem, look at the [original formulation](https://docs.sciml.ai/SciMLSensitivity/stable/examples/pde/pde_constrained/) .
=#

module Example301_PDEOptimSciMLSensitivity
using SkeelBerzins, DifferentialEquations

using DelimitedFiles
using Optimization, OptimizationPolyalgorithms, OptimizationOptimJL, SciMLSensitivity

function f(p)

    ## Problem setup parameters:
    m = 0
    Lx = 10.0
    x = 0.0:0.01:Lx
    dx = x[2] - x[1]

    ## Problem Parameters
    dt = 0.40 * dx^2    # CFL condition
    t0, tMax = 0.0, 1000 * dt
    tspan = (t0, tMax)
    t = t0:dt:tMax

    a0, a1 = p

    function pdefun(x, t, u, dudx)
        c = 1
        f = a1 * dudx
        s = 2.0 * a0 * u

        c, f, s
    end

    icfun(x) = exp(-(x - 3.0)^2)

    function bcfun(xl, ul, xr, ur, t)
        pl = ul
        ql = 0
        pr = ur
        qr = 0

        pl, ql, pr, qr
    end

    params_pdepe = SkeelBerzins.Params(; solver=:DiffEq, nb_design_var=length(p))

    pb = pdepe(m, pdefun, icfun, bcfun, collect(x), tspan; params=params_pdepe)
    prob = DifferentialEquations.ODEProblem(pb)
    sol = DifferentialEquations.solve(prob, ROS34PW1a(; linsolve=SparspakFactorization()); dt=dt, saveat=t)

    sol
end

function main()
    p = [1.0, 1.0] # True solution parameters
    sol_exact = Array(f(p))

    ## Building the Prediction Model
    ps = [0.1, 0.2]  # Initial guess for model parameters
    predict(θ) = Array(f(θ))

    ## Defining Loss function
    function loss(θ)
        pred = predict(θ)
        return sum(abs2.(pred .- sol_exact)), pred # Mean squared error
    end

    LOSS = []     # Loss accumulator
    PRED = []     # prediction accumulator
    PARS = []     # parameters accumulator

    callback = function (θ, l, pred) # callback function to observe training
        display(l)
        append!(PRED, [pred])
        append!(LOSS, l)
        append!(PARS, [θ])
        false
    end

    adtype = Optimization.AutoForwardDiff() # see https://docs.sciml.ai/Optimization/stable/API/optimization_function/#Automatic-Differentiation-Construction-Choice-Recommendations
    optf = Optimization.OptimizationFunction((x, p) -> loss(x), adtype)

    optprob = Optimization.OptimizationProblem(optf, ps)
    res = Optimization.solve(optprob, NewtonTrustRegion(); allow_f_increases=true, callback=callback)

    return res.u
end

using Test

function runtests()
    testval = [1.0, 1.0]
    @test main() ≈ testval
end

end
