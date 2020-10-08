using LinearAlgebra
using Plots
using Revise

includet("IterativeMethods.jl")
includet("InternalArray.jl")
IArray = InternalArray


# Definitions -----------------------------------------------------------------------------

# We take the reference system to be the position of the spring in rest.
# That is, for the one dimensional case: (l, 0).
struct ProblemDefinitions
    t0 # Initial time
    p0 # Initial position
    v0 # Initial speed

    m # Mass
    k # Spring constant
    b # Friction constant
    A # External force amplitude
    ω # External force angular speed
    ProblemDefinitions(; t0 = 0, p0, v0, m, k, b = 0, A = 0, ω = 0) = new(t0, p0, v0, m, k, b, A, ω)
end

function eulermethod_spring(prob::ProblemDefinitions, step, nit)
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

