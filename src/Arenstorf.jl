# Source: https://www.johndcook.com/blog/2020/02/08/arenstorf-orbit/

module Arenstorf

using Plots

mmoon = 0.012277471
mearth = 1-mmoon
p0 = [0.994, 0]
v0 = [0, -2.001585106]

function derivative(t, x)
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

function vderivative(t, x)
    r = zeros(size(x))
    for i in 1:size(x, 2)
        r[:,i] = derivative(t, x[:,i])
    end
    return r
end

function init_plot()
    p = plot(
        cos,
        sin,
        0, 2π,
        label = "moon",
        color = :black,
        style = :dot,
        title = "Arenstorf orbit",
        aspect_ratio = :equal,
        framestyle = :box,
        xlims = (-1.4, 1.2),
        legend = :bottomright,
    )
    scatter!((0, 0), label = "", color = :black, markersize = 1.5)
    return p
end

end

