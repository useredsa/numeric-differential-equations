using Revise

includet("src/Spring.jl")

# 1. Solving the equation using the euler method ------------------------------------------

# We are asked to find solutions for the time interval [t0, t0+40].
# I want to execute the following with a step of 1e-3,
# so I'll set an appropiate number of iterations.
const step = 1e-3
const nit = Int(40/step)  # round to avoid small rounding errors

prob = ProblemDefinitions(p0 = 1.5, v0 = 0, m = 1, b = 0.3, k = 1.5, A = 0.4, ω = 2.4)

x = eulermethod_spring(prob, step, nit)

# 2. Plotting the result ------------------------------------------------------------------

timeline = [prob.t0 + i*step for i in 1:length(x)]
plot(timeline, IArray(x, 1), label = "movement")

savefig("media/spring_std.png")

