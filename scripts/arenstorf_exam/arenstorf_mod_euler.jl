using Revise
using LinearAlgebra
using Plots

includet("src/InternalArray.jl")
includet("src/FixedStepMethods.jl")
includet("src/MixedStepMethods.jl")
includet("src/Interpolation.jl")
includet("src/Arenstorf.jl")
IArray = InternalArray

const closing_err = 5e-3
const stop_time = 40
const mstep = 5e-6

x, time = @timed fixed_step_ode_solver(
    Arenstorf.derivative,
    0,
    [Arenstorf.p0; Arenstorf.v0],
    (t, x) ->
        t > 20 && norm(Arenstorf.p0 .- last(x)[1:2]) < closing_err || t > stop_time,
    step = mstep,
    maxit = 10^7,
    method = mod_euler_step,
)

p = Arenstorf.init_plot()
inds = [i for i in 1:length(x)÷2000:length(x)]
if last(inds) != length(x)
    push!(inds, length(x))
end
plot!(
    IArray(x[inds], 1),
    IArray(x[inds], 2),
    label = "satellite",
    title = "Arenstorf Orbit - Modified Euler's Method",
)
savefig("media/arenstorf/mod_euler.png")

begin
    print("\033[32;1marenstorf>\033[0m \n")
    print("\tMethod: Modified Euler's\n")
    print("\tStep: $(mstep)\n")
    print("\tTarget Error (closing condition): $(closing_err)\n")
    println()
    numit = length(x)-1
    print("\tNumber of iterations: $(numit)\n")
    print("\tRunning time: $(time)s\n")
    print("\tTime (solution): $(numit*mstep)\n")
    print("\tNumber of evaluations: $(2*numit)\n")
    print("\tOrbit closed with target tolerance: $(numit*mstep < stop_time)\n")
    print("\tDistance from last point to start point: ",
          "$(norm(last(x)[1:2] .- Arenstorf.p0))\n")
    n = length(x)
    last5x = [p[1] for p in x[n-4:n]]
    last5y = [p[2] for p in x[n-4:n]]
    closingx = Interpolation.hermite(last5y, last5x, 0)
    print("\tDistance from interpolated closing point to start point: ",
          "$(abs(closingx - Arenstorf.p0[1]))\n")
end

