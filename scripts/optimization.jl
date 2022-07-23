using LinearAlgebra

function bls_steepest(func, grad, α0, x0; rel_error=1e-4, ρ=0.5, c=0.5, return_stats=false)
	α_k = α0
    x_k = x0
    x_prev = zeros(size(x_k, 1))

	typeof(x_k) == Vector{Float64} ? p_k = -Matrix(1.0I, size(x_k, 1), size(x_k, 1)) * grad(x_k) : p_k = -1

	return_stats ? steps = 0 : nothing

    error = 100
	while error > rel_error	
		while func(x_k + α_k.*p_k) > func(x_k) + c*α_k*dot(grad(x_k), p_k)
			α_k = ρ*α_k
		end
        x_prev = x_k
		x_k = x_k + α_k .* p_k

        typeof(x_k) == Vector{Float64} ? error = norm(x_k - x_prev, Inf) / norm(x_k, Inf) : error = abs((x_k - x_prev) / x_k)

		return_stats ? steps = steps + 1 : nothing
	end

	if return_stats == true
		return x_k, error, steps
	else
		return x_k
	end
end

function bls_newton(func, grad, hess, α0, x0; rel_error=1e-4, ρ=0.5, c=0.5, return_stats=false)
	α_k = α0
	x_k = x0
	x_prev = zeros(size(x_k, 1))

	p_k = -hess(x_k) * grad(x_k)
	
	return_stats ? steps = 0 : nothing
	
	error = 100
	while error > rel_error	
		while func(x_k + α_k.*p_k) > func(x_k) + c*α_k*dot(grad(x_k),p_k)
			α_k = ρ*α_k
		end
		x_prev = x_k
		x_k = x_k + α_k .* p_k

		typeof(x_k) == Vector{Float64} ? error = norm(x_k - x_prev, Inf) / norm(x_k, Inf) : error = abs((x_k - x_prev) / x_k)
		
		return_stats ? steps = steps + 1 : nothing
	end
	if return_stats == true
		return x_k, error, steps
	else
		return x_k
	end
end