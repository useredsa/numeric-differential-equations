using Revise

includet("src/Spring.jl")
IArray = InternalArray

# 3. Validating the solution --------------------------------------------------------------

const step = 1e-3
const nit = Int(40/step)  # round to avoid small rounding errors

# We define a problem with no friction
prob = SpringProblem(p0 = 0.8, v0 = 0, m = 1, k = 1.5)

x = eulermethod(prob, step, nit)

timeline = [prob.t0 + i*step for i in 1:length(x)]
plot(timeline, IArray(x, 1), label = "standard")

# 5. Finding the frequency of the unforced movement ---------------------------------------

# Empirically, We can see in the plot that the distance between two
# equal values (!= 0) is aproximately 5.
# Let's get the actual distance from our solution.

p0 = 1       # Reference point
P0 = -1      # Next equal point
for i = 50:length(x) # We avoid the first numbers because they are very close
    if norm(x[p0][1]-x[i][1]) < 1e-3
        P0 = i
        break
    end
end

if P0 != -1
    period = timeline[P0]-timeline[p0] # Is similar to what we observed
    freq = 1/period
end
        

# 6. Using this frequency to force it -----------------------------------------------------

cond_prob = SpringProblem(p0 = 0.8, v0 = 0, m = 1, k = 1.5, A = 1, ω = 2*π*freq)

y = eulermethod(cond_prob, step, nit)
plot!(timeline, IArray(y, 1), label = "conditioned")

savefig("media/spring_frequency_unforced.png")
