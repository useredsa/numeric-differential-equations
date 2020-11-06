using Revise
using LinearAlgebra
using Plots

includet("src/InternalArray.jl")
includet("src/FixedStepMethods.jl")
includet("src/Interpolation.jl")
IArray = InternalArray

G = 6.67430e-11 * 3600^2
earth_mass = 5.9724e-3
sun_mass = 1.9885e3
reduced_mass = 1/(1/sun_mass + 1/earth_mass)
year = 365.242
periheliom = 147.09

r0 = [152.1, 0]
v0 = [0, 2.929e-5] * 3600
step = 150

function f(t, x)
    r = [x[1], x[2]]
    rr = dot(r, r)
    v = [x[3], x[4]]
    return [v; -G*(sun_mass+earth_mass)/rr/sqrt(rr)*r]
end

function full_turn(t, x)
    p2 = x[length(x)-1][1:2]
    p1 = x[length(x)][1:2]
    return p2[2] < 0 && p1[2] >= 0
end

function half_turn(t, x)
    p2 = x[length(x)-1][1:2]
    p1 = x[length(x)][1:2]
    return p2[2] > 0 && p1[2] <= 0
end

methods = [euler_step, mod_euler_step, heum_step, RK4_step]
xs = []
for met in methods
    push!(xs, fixed_step_ode_solver(
        f,
        0,
        [r0; v0],
        full_turn;
        method = met,
        step = step,
        maxit = 10^6
    ))
end

function two_bodies_error(x)
    return maximum(abs.(norm.(x).-norm([r0; v0])))
end

function two_bodies_year(x)
    return (length(x)-1)*step/24
end

plots = []
for i = 1:length(methods)
    print("\033[32;1mtwo_bodies>\033[0m \n")
    print("\tMethod $(methods[i])\n")

    n = length(xs[i])
    last5x = [p[1] for p in xs[i][n-4:n]]
    last5y = [p[2] for p in xs[i][n-4:n]]
    last5t = [Float64(step*i) for i in n-5:n-1]
    closingx = Interpolation.hermite(last5y, last5x, 0)
    closingt = Interpolation.hermite(last5y, last5t, 0)

    m = 2
    while m < n && !half_turn(0, xs[i][1:m])
        m += 1
    end
    half5x = [p[1] for p in xs[i][m-2:m+2]]
    half5y = [p[2] for p in xs[i][m-2:m+2]]
    half5t = [Float64(step*i) for i in m-3:m+1]
    periheliomx = Interpolation.hermite(half5y, half5x, 0)

    print("\tClosing Error: $(norm(last(xs[i])[1:2]-[closingx, 0]))\n")
    print("\tYear: $(closingt/24)\n")
    print("\tYear Error(NASA): $(abs(year-closingt/24))\n")
    print("\tPeriheliom: $(-periheliomx)\n")
    print("\tPeriheliom Error(NASA): $(abs(periheliom+periheliomx))\n\n")

    push!(plots, plot(IArray(xs[i], 1),
                      IArray(xs[i], 2),
                      label = "$(methods[i])",
                      aspect_ratio = :equal))
    savefig("media/two_bodies_$(methods[i]).png")
end

plot(plots[1], plots[2], plots[3], plots[4], title = "Step: $(step)")
savefig("media/two_bodies_step_$(step).png")

