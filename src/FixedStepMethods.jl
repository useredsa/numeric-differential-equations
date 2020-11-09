function euler_step(f, t::Number, x, len::Number)
    return len*f(t, x)
end

function mod_euler_step(f, t::Number, x, len::Number)
    whole_step = len*f(t, x)
    return (whole_step + len*f(t+len, x+whole_step))/2
end

function heum_step(f, t::Number, x, len::Number)
    whole_step = len*f(t, x)
    return (whole_step + 3*len*f(t+2/3*len, x+2/3*whole_step))/4
end

function RK4_step(f, t::Number, x, len::Number)
    k1 = len*f(t, x)
    k2 = len*f(t+len/2, x+k1/2)
    k3 = len*f(t+len/2, x+k2/2)
    k4 = len*f(t+len, x+k3)
    return (k1+2*k2+2*k3+k4)/6
end

macro define_order(f, ord::Integer)
    @eval function $f(s::Symbol)
        if s == :order
            return $ord
        end
        throw(ErrorException("Invalid symbol. Use :order"))
    end
end

@define_order(euler_step, 1)
@define_order(mod_euler_step, 2)
@define_order(heum_step, 2)
@define_order(RK4_step, 4)

function fixed_step_ode_solver(
    f,
    t0::Number,
    x0;
    method=RK4_step,
    step::Number = 1e-3,
    iter::Integer = 10000
)
    t = t0
    x = Array{typeof(x0 + step*f(t0, x0))}(undef, iter+1)
    x[1] = x0
    for j = 1:iter
        x[j+1] = x[j] + method(f, t, x[j], step)
        t += step
    end
    return x
end

function fixed_step_ode_solver(
    f,
    t0::Number,
    x0,
    stop_cond;
    method=RK4_step,
    step::Number = 1e-3,
    maxit::Integer = 10000
)
    t = t0
    x = Array{typeof(x0 + method(f, t, x0, step)), 1}()
    push!(x, x0)
    for j = 1:maxit
        push!(x, last(x) + method(f, t, last(x), step))
        t += step
        if stop_cond(t, x)
            break
        end
    end
    return x
end

