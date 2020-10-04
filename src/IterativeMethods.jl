"""
	eulermethod(f, t0, x0; kwargs...)
	
Solves the differential equation x'(t) = f(t, x(t)) using the euler method.

# Arguments
- `f`: Derivative function
- `t0`: Initial time
- `x0`: Initial position

# Keywords
- `h`: step of the method
- `iter`: number of iterations

# Returns
- `Array{typeof(x0 + h*f(t0, x0))}`: the approximated solution
"""
function eulermethod(f, t0, x0; h = 1e-3, iter = 10000)
	t = t0
	x = Array{typeof(x0 + h*f(t0, x0))}(undef, iter)
	x[1] = x0
	for j = 2:iter
		x[j] = x[j-1] + h*f(t, x[j-1])
		t += h
	end
	return x
end

"""
	targetting_eulermethod(f, dist, t0, x0; kwargs...)
	
Solves the differential equation x'(t) = f(t, x(t))
using a modification of the euler method to stop close to a given target.
The algorithm is similar to the original euler method,
but the function dist(x) specifies the signed distance from the new point to the target,
if it is less than 0 then it is discarded and the method's time step divided by two,
otherwise it becomes part of the solution.
The algorithm stops when the distance is positive and less than the admitted error
or when the maximum iterations limit is reached.

# Arguments
- `f`: Derivative function
- `dist`: Function for the signed distance to the target
- `t0`: Initial time
- `x0`: Initial position

# Keywords
- `h`: step of the method
- `iter`: number of iterations

# Returns
- `Array{typeof(t0)}, Array{typeof(x0 + h*f(t0, x0))}`: times and values obtained
"""
function targetting_eulermethod(f, dist, t0, x0; h = 1e-3, error = 1e-4, maxit = 100000)
	t = typeof(t0+h)[t0]
	x = Array{typeof(x0 + h*f(t0, x0)), 1}()
	push!(x, x0)
	while dist(last(x)) >= error && maxit > 0
		nx = last(x) + h*f(last(t), last(x))
		d = dist(nx)
		if (d < zero(d))
			h /= 2
		else
			push!(t, last(t)+h)
			push!(x, nx)
		end
		maxit -= 1
	end
	return t, x
end

