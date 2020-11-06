module Interpolation

export hermite

function hermite(x, y, x0)
    diff_table = Array{eltype(y), 2}(undef, length(x), length(x))
    diff_table[:,1] = y
    for j = 1:length(x)-1
        for i = 1:length(x)-j
            diff_table[i,j+1] = (diff_table[i+1,j]-diff_table[i,j])/(x[i+j]-x[i])
        end
    end

    y0 = zero(eltype(x0))
    s = 1
    for j = 1:length(x)
        y0 += s*diff_table[1,j]
        s *= x0 .- x[j]
    end
    return y0
end

end # module Interpolation

