using Revise

includet("src/Spring.jl")
IArray = InternalArray


# 3. Validating the solution --------------------------------------------------------------

const step = 1e-3
const nit = Int(40/step)  # round to avoid small rounding errors

# We define a problem with no friction
prob = SpringProblem(p0 = 0.8, v0 = 0, m = 1, k = 1.5)

x = eulermethod(prob, step, nit)

timeline = [prob.t0 + i*step for i in 0:length(x)-1]
plot(timeline, IArray(x, 1), label = "euler method")

# And also define the solution we would have obtained by analytic methods
# Reference: https://en.wikipedia.org/wiki/Simple_harmonic_motion
ω0 = sqrt(prob.k/prob.m)
g(t) = prob.p0*cos(ω0*t) + prob.v0/ω0*sin(ω0*t)

y = g.(timeline)
plot!(timeline, y, label = "analytic solution")

# (The comparison in my computer is very close)
# (Which means that for the second step we could have taken a smaller L very likely)
# Let's check the norm of the difference:

err = maximum(norm.(y-IArray(x, 1)))


# 4. Choosing the step to obtain a final error in the order of 1e-2  ----------------------

# Reference: https://en.wikipedia.org/wiki/Euler_method#Local_truncation_error
# The global truncation error is bounded by $\frac{hM}{2L}(ℯ^{L(t-t_0)}-1)$ where
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

# After checking that with 1e-3 the error is approximately $0.02$
# for the solution without friction.
# It seems natural to admit that an appropiate step to obtain
# an error less than 1e-2 is approximately 5e-4.
# This can be verified using this same script but
# with a smaller error for the method without friction.

# After doing it, I saw that it was not enough and
# the error was still slightly larger than 1e-2.
# With the step equal to 1e-4, the error for the problem without friction is $0.0023$.
# Hence, we can assume that it won't be much larger for the problem with friction.

savefig("media/spring_validation.png")

