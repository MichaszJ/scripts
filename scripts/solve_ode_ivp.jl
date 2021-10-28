function rk_matrices(solver_method)
    if solver_method == "rk_45"
        a = [0, 1/4, 3/8, 12/13, 1, 1/2]
        b = [
            0 0 0 0 0;
            1/4 0 0 0 0;
            3/32 9/32 0 0 0;
            1932/2197 -7200/2197 7296/2197 0 0;
            439/216 -8 3680/513 -845/4104 0;
            -8/27 2 -3544/2565 1859/4104 -11/40
        ]
        c = [16/135, 0, 6656/12825, 28561/56430, -9/50, 2/55]
        cs = [25/216, 0, 1408/2565, 2197/4104, -1/5, 0]
    
    elseif solver_method == "rk_56"
        a = [0, 1/6, 4/15, 2/3, 4/5, 1, 0, 1]
        b = [
            0 0 0 0 0 0 0;
            1/6 0 0 0 0 0 0;
            4/75 16/75 0 0 0 0 0;
            5/6 -8/3 5/2 0 0 0 0;
            -8/5 144/25 -4 16/25 0 0 0;
            361/320 -18/5 407/128 -11/80 55/128 0 0;
            -11/640 0 11/256 -11/160 11/256 0 0;
            93/640 -18/5 803/256 -11/160 99/256 0 1
        ]
        c = [7/1408, 0, 1125/2816, 9/32, 125/768, 0, 5/66, 5/66]
        cs = [31/384, 0, 1125/2816, 9/32, 125/768, 5/66, 0, 0]
    
    elseif solver_method == "rk_78"
        a = [0, 2/27, 1/9, 1/6, 5/12, 1/2, 5/6, 1/6, 2/3, 1/3, 1, 0, 1]
        b = [
            0 0 0 0 0 0 0 0 0 0 0 0;
            2/7 0 0 0 0 0 0 0 0 0 0 0;
            1/36 1/12 0 0 0 0 0 0 0 0 0 0;
            1/24 0 1/82 0 0 0 0 0 0 0 0 0;
            5/12 0 -25/16 25/16 0 0 0 0 0 0 0 0;
            1/20 0 0 1/4 1/5 0 0 0 0 0 0 0;
            -25/108 0 0 125/108 -65/27 125/54 0 0 0 0 0 0;
            31/300 0 0 0 61/225 -2/9 13/900 0 0 0 0 0;
            2 0 0 -53/6 704/45 -107/9 67/90 3 0 0 0 0;
            -91/108 0 0 23/108 -976/135 311/54 -19/60 17/6 -1/12 0 0 0;
            2383/4100 0 0 -341/164 4496/1025 -301/82 2133/4100 45/82 45/164 18/41 0 0;
            3/205 0 0 0 0 -6/41 -3/205 -3/41 3/41 6/41 0 0;
            -1777/4100 0 0 -341/164 4496/1025 -289/82 2193/4100 51/82 33/164 12/41 0 1
        ]
        c = [41/840, 0, 0, 0, 0, 34/105, 9/35, 9/35, 9/280, 9/280, 41/840, 0, 0]
        cs = [0, 0, 0, 0, 0, 34/105, 9/35, 9/35, 9/280, 9/280,0, 41/840, 41/480]
    end
    
    return a, b, c, cs
end

function rk_solver(ode_func, y_init, t_init, t_final; solver_method="rk_45", tolerance=100.0, beta=0.8)
    a, b, c, cs = rk_matrices(solver_method)

    step_size = (t_final - t_init) / 100

    t_out = convert(Array{Float64}, [t_init])
    
    y_size = size(y_init, 1)
    y_out = [[y_init[i]] for i=1:size(y_init, 1)]
    
    f_matrix = ones(size(a, 1), size(y_init, 1))

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

            step_size = step_size * beta * (tolerance / truncation_err) ^ (1/(size(a, 1) - 1))
        end        
    end
    return y_out, t_out
end