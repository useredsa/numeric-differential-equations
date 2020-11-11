using LinearAlgebra

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
    method = mod_euler_step,
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
        aux = last(x) .+ method(f, last(t), last(x), [step, step/2])
        xfull = aux[:,1]
        xhalf = aux[:,2] + method(f, last(t)+step/2, aux[:,2], step/2)

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
        step = min(q*step, max_step)
        if step < min_step
            throw(ErrorException("step went below the min_step($(min_step)) value"))
        end
    end
    return t, x
end

struct ExtButcherTableau
    a
    c
    b1
    b2
    order

    ExtButcherTableau(; a, c, b1, b2, order::Integer) = new(a, c, b1, b2, order)
end

DORPRI5_tableau = ExtButcherTableau(
    a = [          0  1/5        3/40    44/45     19372/6561   9017/3168     35/384;
                   0    0        9/40   -56/15    -25360/2187     -355/33          0;
                   0    0           0     32/9     64448/6561  46732/5247   500/1113;
                   0    0           0        0       -212/729      49/176    125/192;
                   0    0           0        0              0 -5103/18656 -2187/6784;
                   0    0           0        0              0           0      11/84;
                   0    0           0        0              0           0          0],
    c = [          0, 1/5,       3/10,     4/5,           8/9,          1,         1],

    b1 = [    35/384,   0,   500/1113, 125/192,    -2187/6784,      11/84,         0],
    b2 = [5179/57600,   0, 7571/16695, 393/640, -92097/339200,   187/2100,      1/40],
    order = 5
)

RKFehlberg_tableau = ExtButcherTableau(
    a = [       0  1/4        3/32    1932/2197   439/216      -8/27;
                0    0        9/32   -7200/2197        -8          2;
                0    0           0    7296/2197  3680/513 -3544/2565;
                0    0           0            0 -845/4104  1859/4104;
                0    0           0            0         0     -11/40;
                0    0           0            0         0          0],
    c = [       0, 1/4,        3/8,       12/13,        1,       1/2],

    b1 = [ 16/135,   0, 6656/12825, 28561/56430,    -9/50,      2/55],
    b2 = [ 25/216,   0,  1408/2565,   2197/4104,     -1/5,         0],
    order = 4
)

function double_step(tab::ExtButcherTableau, f, t, x, step)
    k = Array{typeof(x + step*f(t, x)), 1}(undef, length(tab.c))
    k[1] = step*f(t + step*tab.c[1], x)
    for j = 2:length(k)
        k[j] = step*f(t + step*tab.c[j], x + tab.a[1:(j-1), j]' * k[1:(j-1)])
    end
    return (x + tab.b1' * k, x + tab.b2' * k)
end

# function double_step(

function adaptative_rk_ode_solver(
    f,
    t0::Number,
    x0,
    stop_cond;
    tableau::ExtButcherTableau = RKFehlberg_tableau,
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
        x1, x2 = double_step(tableau, f, last(t), last(x), step)

        err = norm(x1 .- x2)
        if err < target_err
            push!(t, last(t)+step)
            push!(x, x1)
            if stop_cond(t, x)
                break
            end
        end

        q = (target_err*step/err) ^ (1/(tableau.order))
        q = min(4, max(0.1, q))
        step = min(q*step, max_step)
        if step < min_step
            throw(ErrorException("step went below the min_step($(min_step)) value"))
        end
    end
    return t, x
end

