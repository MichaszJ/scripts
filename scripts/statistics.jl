function mean(x)
    return sum(x) / size(x, 1)
end

function trimmed_mean(x, p)
    return sum(x[p + 1 : end - p]) / (size(x, 1) - 2*p)
end

function weighted_mean(x, w)
    return sum(x .* w) / sum(w)
end

function median(x)
    x_sorted = sort(x)
    n = size(x_sorted,1)
    return n % 2 == 0 ? (x_sorted[n ÷ 2] + x_sorted[(n ÷ 2) + 1]) / 2 : x_sorted[(n+1) ÷ 2]
end

function weighted_median(x, w)
    perm = sortperm(x)
    x_sorted = x[perm]
    w_sorted = w[perm]

    n = size(x_sorted, 1)
    cumulative_weights = [sum(w_sorted[1:i]) for i in 1:n]

    comp_vec = [cumulative_weights[i] > 0.5*sum(w_sorted) for i in 1:n]

    index = argmax(comp_vec)

    weighted_vec = x_sorted .* w_sorted
    sum(weighted_vec[1:index]) == 0.5*sum(weighted_vec) ? weighted_median = x_sorted[index] : weighted_median = (x_sorted[index] + x_sorted[index+1])/2

    return weighted_median
end

function mean_absolute_deviation(x)
    x̄ = mean(x)
    return sum((x .- x̄).^2) / size(x, 1)
end

function variance(x)
    x̄ = mean(x)
    return sum((x .- x̄).^2) / (size(x, 1) - 1)
end

function standard_deviation(x)
    return √variance(x)
end

function median_absolute_deviation(x)
    med = median(x)
    return median(abs.(x .- med))
end

function quartiles(x)
    n = size(x, 1)
    x_sorted = sort(x)
    
    first_quart = median(x_sorted[1:n ÷ 2])
    second_quart = median(x_sorted)
    
    if n % 2 == 0
        third_quart = median(x_sorted[n÷2 + 1: n])
    else
        third_quart = median(x_sorted[n÷2 + 2 : n])
    end
        
    return convert.(Int64, [first_quart, second_quart, third_quart])
end

function interquartile_range(x)
    q = quartiles(x)
    return q[3] - q[1]
end