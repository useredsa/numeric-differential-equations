# Source: https://www.johndcook.com/blog/2020/02/08/arenstorf-orbit/

using Revise
using LinearAlgebra
using Plots

includet("src/InternalArray.jl")
includet("src/FixedStepMethods.jl")
includet("src/MixedStepMethods.jl")
IArray = InternalArray

mmoon = 0.012277471
mearth = 1-mmoon
p0 = [0.994, 0]
v0 = [0, -2.001585106]
target_err = 1e-4

function derivative(t, x)
    if ndims(x) == 2
        r = zeros(size(x))
        for i in 1:size(x, 2)
            r[:,i] = derivative(t, x[:,i])
        end
        return r
    end
    D1 = ((x[1]+mmoon)*(x[1]+mmoon) + x[2]*x[2])
    D1 *= sqrt(D1)
    D2 = ((x[1]-mearth)*(x[1]-mearth) + x[2]*x[2])
    D2 *= sqrt(D2)
    return [
        x[3],
        x[4],
        x[1] + 2*x[4] - mearth*(x[1]+mmoon)/D1 - mmoon*(x[1]-mearth)/D2,
        x[2] - 2*x[3] - mearth*x[2]/D1 - mmoon*x[2]/D2
    ]
end

t, x = mixed_step_ode_solver(
    derivative,
    0,
    [p0; v0],
    (t, x) ->
        last(t) > 5 && norm(p0 .- last(x)[1:2]) < target_err*last(t) || last(t) > 50,
    min_step = 1e-7,
    max_step = 0.1,
    target_err = target_err,
    method = mod_euler_step,
);

plot(
    IArray(x, 1),
    IArray(x, 2),
    label = "satellite",
    title = "Arenstorf orbit",
    aspect_ratio = :equal,
    framestyle = :box,
    xlims = (-1.4, 1.2),
    legend = :bottomright
)
plot!(cos, sin, 0, 2π, label = "moon", color = :black, style = :dot)
scatter!((0, 0), label = "", color = :black, markersize = 1.5)
# savefig("media/arenstorf_orbit.png")

using IJulia

px = [x[i][1] for i in 1:length(x)]
py = [x[i][2] for i in 1:length(x)]
anim = Animation()
i = 1
for j = 0:0.1:last(t)
    while i < length(t) && t[i] < j
        i += 1
    end
    plot(
        px[1:i],
        py[1:i],
        label = "satellite",
        title = "Arenstorf orbit",
        aspect_ratio = :equal,
        framestyle = :box,
        xlims = (-1.4, 1.2),
        ylims = (-1.3, 1.4),
        legend = :bottomright,
    )
    inds = range(0, 2*π*j/last(t), length = 150)
    plot!(cos.(inds), sin.(inds), label = "moon", color = :black, style = :dot)
    scatter!((0, 0), label = "", color = :black, markersize = 1.5)
    frame(anim)
end
gif(anim, "media/arenstorf_orbit.gif", fps = 15)

