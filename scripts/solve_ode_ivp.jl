function rk_45(ode_func, y_init, t_init, t_final, tolerance=0.2, beta=0.8)
    a = [0, 1/4, 3/8, 12/13, 1, 1/2]
    b = [
        0 0 0 0 0;
        1/4 0 0 0 0;
        3/32 9/32 0 0 0;
        1932/2197 -7200/2197 7296/2197 0 0;
        439/216 -8 3680/513 -845/4104 0;
        -8/27 2 -3544/2565 1859/4104 -11/40
    ]
    cs = [25/216, 0, 1408/2565, 2197/4104, -1/5, 0]
    c = [16/135, 0, 6656/12825, 28561/56430, -9/50, 2/55]

    step_size = (t_final - t_init) / 100

    t_out = convert(Array{Float64}, [t_init])
    
    y_size = size(y_init, 1)
    y_out = [[y_init[i]] for i=1:size(y_init, 1)]
    
    f_matrix = ones(6, size(y_init, 1))

    while t_init < t_final
        for i = 1:size(f_matrix, 1)
            t_delta = a[i] * step_size
            y_delta = step_size .* broadcast(sum, [b[[i], :] .* f_matrix[:5, j] for j=1:y_size])
            
            f_matrix[[i],:] = ode_func(t_init + t_delta, y_init .+ y_delta)
        end

        truncation_vec = broadcast(abs, step_size * (c - cs) .* f_matrix)
        truncation_err = maximum(truncation_vec)

        if truncation_err > tolerance
            step_size = step_size * beta * (tolerance / truncation_err)^(1/5)
        else
            y_init = y_init .+ step_size * broadcast(sum, [c .* f_matrix[:, j] for j=1:y_size])    
            t_init = t_init .+ step_size

            for i=1:size(y_init, 1)
                push!(y_out[i], y_init[i])
            end

            push!(t_out, t_init)

            step_size = step_size * beta * (tolerance / truncation_err) ^ (1/5)
        end        
    end
    return y_out, t_out
end