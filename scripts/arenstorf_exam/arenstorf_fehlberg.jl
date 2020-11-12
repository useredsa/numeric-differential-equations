using Revise
using LinearAlgebra
using Plots

includet("src/InternalArray.jl")
includet("src/FixedStepMethods.jl")
includet("src/MixedStepMethods.jl")
includet("src/Interpolation.jl")
includet("src/Arenstorf.jl")
IArray = InternalArray

const target_err = 1e-12
const closing_err = 1e-5
const stop_time = 40
const min_step = 1e-7
const max_step = 0.1
extra_info = Dict()

(t, x), time = @timed adaptative_rk_ode_solver(
    Arenstorf.derivative,
    0,
    [Arenstorf.p0; Arenstorf.v0],
    (t, x) ->
        last(t) > 20 && norm(Arenstorf.p0 .- last(x)[1:2]) < closing_err || last(t) > stop_time,
    min_step = min_step,
    max_step = max_step,
    target_err = target_err,
    tableau = RKFehlberg_tableau,
    maxit = 5*10^6,
    debug_dict = extra_info
)

steps = [t[i]-t[i-1] for i in 2:length(t)]

p1 = Arenstorf.init_plot()
plot!(
    IArray(x, 1),
    IArray(x, 2),
    label = "satellite",
    title = "Arenstorf Orbit - RK-Fehlberg",
)
p2 = plot(
    t[2:length(t)],
    steps,
    label = "Step Size",
    title = "Variation of step",
    legend = :bottomright,
)
plot(p1, p2)
savefig("media/arenstorf/fehlberg.png")

begin
    print("\033[32;1marenstorf>\033[0m \n")
    print("\tMethod: Felhberg's\n")
    print("\tMin. Step: $(min_step)\n")
    print("\tMax. Step: $(max_step)\n")
    print("\tTolerance: $(target_err)\n")
    print("\tTarget Error (closing condition): $(closing_err)\n")
    println()
    numit = extra_info[:numit]
    print("\tNumber of iterations: $(numit)\n")
    print("\tRunning time: $(time)s\n")
    print("\tTime (solution): $(last(t))\n")
    print("\tNumber of evaluations: $(6*numit)\n")
    print("\tOrbit closed with target tolerance: $(last(t) < stop_time)\n")
    print("\tDistance from last point to start point: ",
          "$(norm(last(x)[1:2] .- Arenstorf.p0))\n")
    n = length(x)
    last5x = [p[1] for p in x[n-4:n]]
    last5y = [p[2] for p in x[n-4:n]]
    closingx = Interpolation.hermite(last5y, last5x, 0)
    print("\tDistance from interpolated closing point to start point: ",
          "$(abs(closingx - Arenstorf.p0[1]))\n")
end


# Creation of Gif

using IJulia

len = length(t)÷2
p1 = plot(
    title = "Arenstorf orbit",
    aspect_ratio = :equal,
    framestyle = :box,
    xlims = (-1.4, 1.2),
    ylims = (-1.2, 1.2),
    legend = :bottomright,
)
plot!(p1, [x[1][1]], [x[1][2]], label = "satellite", style = :auto)
plot!(p1, [cos(0)], [sin(0)], color = :black, label = "moon")
scatter!((0, 0), label = "", color = :black, markersize = 1.5)
p2 = plot(
    title = "Variation of step",
    xlims = (0, t[len]),
    ylims = (0, maximum(steps[1:len-1])*1.2),
    legend = :bottomright,
)
plot!(p2, [0], [0], label = "step size", style = :auto)


anim = Animation()
i = 2
for j = 0:0.05:t[len]
    while i < len && t[i] < j
        i += 1
    end
    push!(p1, [(x[i][1], x[i][2]), (cos(2*π*t[i]/t[len]), sin(2*π*t[i]/t[len])), (NaN, NaN)])
    push!(p2, t[i], steps[i-1])
    plot(p1, p2)
    frame(anim)
end
gif(anim, "media/arenstorf/fehlberg.gif", fps = 15)

