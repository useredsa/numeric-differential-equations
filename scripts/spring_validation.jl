using Revise

includet("src/Spring.jl")

# 3. Validating the solution --------------------------------------------------------------

const step = 1e-3
const nit = Int(40/step)  # round to avoid small rounding errors

# We define a problem with no friction
prob = ProblemDefinitions(p0 = 1.5, v0 = 0, m = 1, k = 1.5)

x = eulermethod_spring(prob, step, nit)

timeline = [prob.t0 + i*step for i in 1:length(x)]
plot(timeline, IArray(x, 1), label = "euler method")

# And also define the solution we would have obtained by analytic methods
ω0 = sqrt(prob.k/prob.m)
g(t) = prob.p0*cos(ω0*t) + prob.v0/ω0*sin(ω0*t)

y = g.(timeline)
plot!(timeline, y, label = "analytic solution")

# (The comparison in my computer is very close)
# (Which means that for the second step we could have taken a smaller L very likely)


# 4. Choosing the step to obtain a final error in the order of 10e2  ----------------------

# Reference: https://en.wikipedia.org/wiki/Euler_method#Local_truncation_error
# The global truncation error is bounded by $\frac{hM}{2L}(ℯ^(L(t-t_0))-1)$ where
# M is an upper bound for the second derivative of the solution (in the interval) and
# L is the lipschitz constant of f.

# This will be explained later, but L < k+b and M can be taken to be maximum of the
# second derivative for the approximated solution.

# <-- Copied from the definitions
n = 1
function f(t, x)
    p = x[1:n]
    v = x[n+1:2*n]
    return [v; (-prob.k*p -prob.b*v + [prob.A*sin(prob.ω*t)])/prob.m]
end

M = maximum((norm ∘ f).(timeline, x))
L = prob.k+prob.b
err_ct = M/L*(ℯ^(L*40)-1) # error constant
target_step = 10e-2/err_ct

savefig("media/spring_validation.png")

