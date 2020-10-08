using Revise

includet("src/Spring.jl")
IArray = InternalArray


# 1. Solving the equation using the euler method ------------------------------------------

# We are asked to find solutions for the time interval [t0, t0+40].
# I want to execute the following with a step of 1e-3,
# so I'll set an appropiate number of iterations.
const step = 1e-3
const nit = Int(40/step)  # round to avoid small rounding errors

prob = SpringProblem(p0 = 0.8, v0 = 0, m = 1, k = 1.5, b = 0.3, A = 0.4, ω = 2.4)

x = eulermethod(prob, step, nit)

# 2. Plotting the result ------------------------------------------------------------------

timeline = [prob.t0 + i*step for i in 0:length(x)-1]
position = plot(timeline, IArray(x, 1), label = "movement")
speed = plot(timeline, IArray(x, 2), label = "speed")
plot(position, speed)

savefig("media/spring_std.png")

