# Example 301: PDE Constrained Optimization

This example is used to show the integration of the package in the SciML ecosystem.
For more details on the problem, look at the following link: [https://docs.sciml.ai/SciMLSensitivity/stable/examples/pde/pde_constrained/](https://docs.sciml.ai/SciMLSensitivity/stable/examples/pde/pde_constrained/).


```
module Example301_SciMLSensitivityIntegration
using SkeelBerzins

using DelimitedFiles
using Optimization, OptimizationPolyalgorithms, OptimizationOptimJL, SciMLSensitivity

function f(p)

    # Problem setup parameters:
    m = 0
    Lx = 10.0
    x  = 0.0:0.01:Lx
    dx = x[2] - x[1]

    # Problem Parameters
    dt       = 0.40*dx^2    # CFL condition
    t0, tMax = 0.0 ,1000*dt
    tspan    = (t0,tMax)
    t        = t0:dt:tMax

    a0, a1 = p

    function pdefun(x,t,u,dudx)
        c = 1
        f = a1*dudx
        s = 2.0*a0*u

        c,f,s
    end

    icfun(x) = exp(-(x-3.0)^2)

    function bcfun(xl,ul,xr,ur,t)
        pl = ul 
        ql = 0
        pr = ur
        qr = 0

        pl,ql,pr,qr
    end

    params_pdepe = SkeelBerzins.Params()
    params_pdepe.solver = :DiffEq

    pb = pdepe(m,pdefun,icfun,bcfun,collect(x),tspan ; params=params_pdepe)
    prob = DifferentialEquations.ODEProblem(pb)
    sol = DifferentialEquations.solve(prob,RadauIIA3(linsolve = SparspakFactorization()), dt=dt,saveat=t)

    sol
end

function main()
    p       = [1.0,1.0] # True solution parameters
    sol_exact = f(p)

    # Building the Prediction Model
    ps  = [0.1, 0.2]  # Initial guess for model parameters
    function predict(θ)
        sol = f(θ)
        Array(sol)
    end

    # Defining Loss function
    function loss(θ)
        pred = predict(θ)
        l = predict(θ)  - sol_exact
        return sum(abs2, l), pred # Mean squared error
    end

    LOSS  = []     # Loss accumulator
    PRED  = []     # prediction accumulator
    PARS  = []     # parameters accumulator

    callback = function (θ,l,pred) # callback function to observe training
        # display(l)
        append!(PRED, [pred])
        append!(LOSS, l)
        append!(PARS, [θ])
        false
    end


    adtype = Optimization.AutoForwardDiff() # see https://docs.sciml.ai/Optimization/stable/API/optimization_function/#Automatic-Differentiation-Construction-Choice-Recommendations
    optf = Optimization.OptimizationFunction((x,p)->loss(x), adtype)

    optprob = Optimization.OptimizationProblem(optf, ps)
    res = Optimization.solve(optprob, NewtonTrustRegion(), allow_f_increases=true, callback=callback)

    return res.u
end

function test()
    testval = [1.0,1.0]
    main() ≈ testval
end

end
```
