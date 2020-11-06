using LinearAlgebra
using Plots
using Revise

includet("src/IterativeMethods.jl")
includet("src/InternalArray.jl")

IArray = InternalArray

t0 = 0
x0 = [0, 100]
v0 = [300, 0]
m = 1
c = 0.1
g = [0, -9.8]

# Solving without friction
x_nofriction = eulermethod(t0, x0) do t, x
	v0 + t*g
end
plot_nofriction = plot(IArray(x_nofriction, 1), IArray(x_nofriction, 2), label = "without friction")

# Solving with friction
x_friction = eulermethod(t0, [x0; v0]) do t, x
	v = [x[3], x[4]]
	return [v; g-c/m*norm(v)*v]
end
plot_friction = plot(IArray(x_friction, 1), IArray(x_friction, 2), label = "with friction")

plot(plot_nofriction, plot_friction, plot_title = "comparing with and without friction")


# Solving for collision with floor
t, x = targetting_eulermethod(x -> x[2], t0, [x0; v0]; h = 5e-3) do t, x
	v = [x[3], x[4]]
	return [v; g-c/m*norm(v)*v]
end
p = plot(IArray(x, 1), IArray(x, 2), label = "landing")



# Needed to stop showing due to error
using IJulia

x1 = IArray(x, 1)
x2 = IArray(x, 2)
v1 = IArray(x, 3)
v2 = IArray(x, 4)

ppp = 35 # points per frame
lims(v) = (minimum(v), maximum(v))
p1 = plot(1, xlim = lims(x1), ylim = lims(x2), label = "position")
p2 = plot(1, xlim = lims(v1), ylim = lims(v2), label = "speed")
anim = Animation()
for i = ppp:ppp:length(x)
	for j = i-ppp+1:i
    	push!(p1, x1[j], x2[j])
    	push!(p2, v1[j], v2[j])
    end
	plot(p1, p2)
	frame(anim)
end
gif(anim, "media/dropped-object-trayectory-to-collision.gif", fps = 24)

