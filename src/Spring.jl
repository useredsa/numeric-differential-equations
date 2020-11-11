using LinearAlgebra
using Plots
using Revise

includet("IterativeMethods.jl")
includet("InternalArray.jl")

# Definitions -----------------------------------------------------------------------------

# We take the reference system to be the position of the spring in rest.
# That is, for the one dimensional case: (l, 0).

"""
    SpringProblem
    
Contains the parameters associated with a typical spring problem.

# Fields
- `t0`: Initial time of the problem (defaults to 0)
- `p0`: Initial position
- `x0`: Initial speed

- `m`: Mass of the object
- `k`: Spring's constant
- `b`: Friction constant
- `A`: External force amplitude
- `ω`: External force angular speed
"""
struct SpringProblem
    t0
    p0
    v0

    m
    k
    b
    A
    ω
    
    SpringProblem(; t0 = 0, p0, v0, m, k, b = 0, A = 0, ω = 0) = new(t0, p0, v0, m, k, b, A, ω)
end

function eulermethod(prob::SpringProblem, step, nit)
    # We ought to transform the equation to a  differential equation of order one,
    # for which we write x = [position; speed].
    # That is, the variable x will represent the solution rather than the position.
    # And the solution contains the position and the speed.
    n = length(prob.p0)
    function f(t, x)
        p = x[1:n]
        v = x[n+1:2*n]
        return [v; (-prob.k*p -prob.b*v + [prob.A*sin(prob.ω*t)])/prob.m]
    end
    

    return eulermethod(f, prob.t0, [prob.p0; prob.v0], h = step, iter = nit)
end

function eulermethod_by_approximation(prob::SpringProblem, timelen; error = 1e-3)
    n = length(prob.p0)
    function f(t, x)
        p = x[1:n]
        v = x[n+1:2*n]
        return [v; (-prob.k*p -prob.b*v + [prob.A*sin(prob.ω*t)])/prob.m]
    end
    

    return eulermethod_by_approximation(f, prob.t0, timelen, [prob.p0; prob.v0]; error = error)
end

