using Revise
using LinearAlgebra
using Plots

includet("src/InternalArray.jl")
includet("src/FixedStepMethods.jl")
includet("src/MixedStepMethods.jl")
includet("src/Arenstorf.jl")
IArray = InternalArray

target_err = 1e-4
stop_time = 20
mstep = 1e-5

x = fixed_step_ode_solver(
    Arenstorf.derivative,
    0,
    [Arenstorf.p0; Arenstorf.v0],
    (t, x) ->
        t > 5 && norm(Arenstorf.p0 .- last(x)[1:2]) < target_err*t || t > stop_time,
    step = mstep,
    maxit = 5*10^6,
    method = euler_step,
)

p = Arenstorf.init_plot()
inds = 1:length(x)÷200:length(x)
plot!(
    IArray(x[inds], 1),
    IArray(x[inds], 2),
    label = "satellite",
    title = "Arenstorf Orbit - Euler's Method",
)
savefig("media/arenstorf/euler.png")

begin
    print("\033[32;1mtwo_bodies>\033[0m \n")
    print("\tMethod: Euler's\n")
    print("\tStep: $(mstep)\n")
    print("\tTarget Error (closing condition): $(target_err)\n")
    println()
    numit = length(x)-1
    print("\tNumber of iterations: $(numit)\n")
    print("\tTime (solution): $(numit*mstep)\n")
    print("\tNumber of evaluations: $(numit)\n")
    print("\tOrbit closed with target tolerance: $(numit*mstep < stop_time)\n")
    print("\tDistance from last point to start point: ",
          "$(norm(norm(last(x)[1:2] .- Arenstorf.p0)/(mstep*numit)))\n")
end

