function err_estimation(xfull, xhalf, order::Integer)
    return (1<<order)/(1<<order-1)*norm(xhalf .- xfull)
end

function richardson_extrapolation(xfull, xhalf, order::Integer)
    return ((1<<order)*xhalf .- xfull)/(1<<order-1)
end

function mixed_step_ode_solver(
    f,
    t0::Number,
    x0,
    stop_cond;
    method=mod_euler_step,
    min_step::AbstractFloat = 1e-3,
    max_step::AbstractFloat = 10,
    target_err::AbstractFloat = 0.1,
    step0::AbstractFloat = 1.0,
    maxit::Integer = 10000
)
    t = Array{typeof(t0 + step0), 1}()
    x = Array{typeof(x0 + step0*f(t0, x0)), 1}()
    step = step0
    push!(t, t0)
    push!(x, x0)
    i = 1
    for j = 1:maxit
        xfull = last(x) + method(f, last(t), last(x), step)
        xhalf = last(x) + method(f, last(t), last(x), step/2)
        xhalf = xhalf + method(f, last(t)+step/2, xhalf, step/2)

        err_full = err_estimation(xfull, xhalf, method(:order))
        if err_full < target_err
            push!(t, last(t)+step)
            push!(x, richardson_extrapolation(xfull, xhalf, method(:order)))
            if stop_cond(t, x)
                break
            end
        end

        q = (target_err*step/err_full) ^ (1/(method(:order)))
        q = min(4, max(0.1, q))
        # println(xfull, xhalf)
        # print(last(t), " ", last(x), " ", err_full, " ", q, "\n")
        step = min(q*step, max_step)
        if step < min_step
            throw(ErrorException("step went below the min_step($(min_step)) value"))
        end
    end
    return t, x
end

