using Revise

includet("src/Spring.jl")
IArray = InternalArray


# 1. Solving the equation using the euler method by approximation ------------------------

const timelen = 40

prob = SpringProblem(p0 = 0.8, v0 = 0, m = 1, k = 1.5, b = 0.3, A = 0.4, ω = 2.4)

x = eulermethod_by_approximation(prob, timelen, error=1e-4)

# 2. Plotting the result ------------------------------------------------------------------

timeline = [prob.t0+timelen*i/length(x) for i in 0:length(x)-1]
pos = plot(timeline, IArray(x, 1), label = "movement")
speed = plot(timeline, IArray(x, 2), label = "speed")
plot(pos, speed)

savefig("media/spring_by_approximation.png")

