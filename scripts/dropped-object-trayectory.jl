using LinearAlgebra
using Plots
using Revise

includet("src/IterativeMethods.jl")
includet("src/InternalArray.jl")

x0 = [0, 300]
v0 = [100, 0]
m = 1
c = 0.1
g = [0, -9.8]

# Solving for speed
v = eulermethod(0, v0) do t, v
	g
end
p1 = plot(InternalArray(v, 1), InternalArray(v, 2), label = "v (no)")

# Solving por position
x = eulermethod(0, x0) do t, x
	v0 + t*g
end
p2 = plot(InternalArray(x, 1), InternalArray(x, 2), label = "x (no)")

# Solving for speed with friction
v = eulermethod(0, v0) do t, v
	g-c/m*norm(v)*v
end
p3 = plot(InternalArray(v, 1), InternalArray(v, 2), label = "v (yes)")

# Solving for position with friction
x = eulermethod(0, [x0; v0]) do t, x
	v = [x[3], x[4]]
	return [v; g-c/m*norm(v)*v]
end
p4 = plot(InternalArray(x, 1), InternalArray(x, 2), label = "x (yes)")

plot(p1, p2, p3, p4, plot_title = "Position and speed with and without friction")

