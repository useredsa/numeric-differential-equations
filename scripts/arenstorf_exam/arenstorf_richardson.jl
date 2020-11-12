using Revise
using LinearAlgebra
using Plots

includet("src/InternalArray.jl")
includet("src/FixedStepMethods.jl")
includet("src/MixedStepMethods.jl")
includet("src/Arenstorf.jl")
IArray = InternalArray

target_err = 1e-10
closing_err = 1e-5
stop_time = 40
min_step = 1e-7
max_step = 0.1
extra_info = Dict()

t, x = mixed_step_ode_solver(
    Arenstorf.vderivative,
    0,
    [Arenstorf.p0; Arenstorf.v0],
    (t, x) ->
        last(t) > 20 && norm(Arenstorf.p0 .- last(x)[1:2]) < closing_err || last(t) > stop_time,
    min_step = min_step,
    max_step = max_step,
    target_err = target_err,
    method = RK4_step,
    maxit = 5*10^6,
    debug_dict = extra_info
);

steps = [t[i]-t[i-1] for i in 2:length(t)]

p1 = Arenstorf.init_plot()
plot!(
    IArray(x, 1),
    IArray(x, 2),
    label = "satellite",
    title = "Arenstorf Orbit - RK4-Richardson",
)
p2 = plot(
    t[2:length(t)],
    steps,
    label = "Step Size",
    title = "Variation of step",
    legend = :bottomright,
)
plot(p1, p2)
savefig("media/arenstorf/richardson.png")

begin
    print("\033[32;1mtwo_bodies>\033[0m \n")
    print("\tMethod: Richardson of RK4\n")
    print("\tMin. Step: $(min_step)\n")
    print("\tMax. Step: $(max_step)\n")
    print("\tTarget Error (closing condition): $(target_err)\n")
    println()
    numit = extra_info[:numit]
    print("\tNumber of iterations: $(numit)\n")
    print("\tTime (solution): $(last(t))\n")
    print("\tNumber of evaluations: $(11*numit)\n")
    print("\tOrbit closed with target tolerance: $(last(t) < stop_time)\n")
    print("\tDistance from last point to start point: ",
          "$(norm(norm(last(x)[1:2] .- Arenstorf.p0)/last(t)))\n")
end

# Creation of Gif

using IJulia

px = [x[i][1] for i in 1:length(x)]
py = [x[i][2] for i in 1:length(x)]
p1 = plot(
    label = ["moon", "satellite"],
    color = :black,
    style = [:dot, :dash],
    title = "Arenstorf orbit",
    aspect_ratio = :equal,
    framestyle = :box,
    xlims = (-1.4, 1.2),
    ylims = (-1.2, 1.2),
    legend = :bottomright,
)
scatter!((0, 0), label = "", color = :black, markersize = 1.5)
plot!(p1, [px[1]], [py[1]], style = :dash, label = "satellite")
plot!(p1, [cos(0)], [sin(0)], style = :dot, label = "moon")
p2 = plot(
    title = "Variation of step",
    xlims = (0, last(t)),
    ylims = (0, 0.012),
    legend = :bottomright,
)
plot!(p2, [0], [1], label = "step size")
    
    
anim = Animation()
i = 2
len = length(t)÷2
for j = 0:0.05:t[len]
    while i < len && t[i] < j
        i += 1
    end
    push!(p1, [(cos(2*π*t[i]/t[len]), sin(2*π*t[i]/t[len])), (px[i], py[i])])
    push!(p2, t[i], steps[i-1])
    plot(p1, p2)
    frame(anim)
end
gif(anim, "media/arenstorf/richardson.gif", fps = 15)

