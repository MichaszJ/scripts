using LinearAlgebra

function get_orbital_elements(r⃗, v⃗; μ=398600, elements_type="standard")
    r = norm(r⃗)
    v = norm(v⃗)

    h⃗ = r⃗ × v⃗
    h = norm(h⃗)

    i = acos(h⃗[3] / h)

    N⃗ = [0, 0, 1] × h⃗
    N = norm(N⃗)

    N⃗[2] >= 0 ? Ω = acos(N⃗[1] / N) : Ω = 2 * pi - acos(N⃗[1] / N)

    v_rad = (r⃗ ⋅ v⃗) / r
    e⃗ = (1 / μ) * ((v^2 - μ / r) * r⃗ - r * v_rad * v⃗)
    e = norm(e⃗)

    e⃗[3] >= 0 ? ω = acos(dot(N⃗ / N, e⃗ / e)) : ω = 2 * π - acos(dot(N⃗ / N, e⃗ / e))

    v_rad >= 0 ? θ = acos(dot(e⃗ / e, r⃗ / r)) : θ = 2 * π - acos(dot(e⃗ / e, r⃗ / r))

    a = (r * (1 + e * cos(θ))) / (1 - e^2)

    if elements_type == "standard"
        return [a, e, i, Ω, ω, θ]
    elseif elements_type == "curtis"
        return [h, i, Ω, e, ω, θ]
    end
end

function gibbs_method(r⃗₁, r⃗₂, r⃗₃; μ=398600, error_tolerance=1e-5, return_state_vector=false, elements_type="standard")
	r₁, r₂, r₃ = norm(r⃗₁), norm(r⃗₂), norm(r⃗₃)
	C⃗₁₂, C⃗₂₃, C⃗₃₁ = r⃗₁ × r⃗₂, r⃗₂ × r⃗₃, r⃗₃ × r⃗₁

    coplanar_error = abs((r⃗₁ ./ r₁) ⋅ (C⃗₂₃ ./ norm(C⃗₂₃)))
	if coplanar_error > error_tolerance
		error("Error: Vectors are not coplanar, |ûᵣ₁ ⋅ Ĉ₂₃| = $(round(coplanar_error, digits=8)) > $(error_tolerance)\nAdjust error_tolerance or choose different position vectors")
	
	else
		N⃗ = r₁*C⃗₂₃ + r₂*C⃗₃₁ + r₃*C⃗₁₂
		D⃗ = C⃗₁₂ + C⃗₂₃ + C⃗₃₁
		S⃗ = r⃗₁.*(r₂ - r₃) + r⃗₂.*(r₃ - r₁) + r⃗₃.*(r₁ - r₂)

		v⃗₂ = sqrt(μ / (norm(N⃗) * norm(D⃗))) * ((D⃗ × r⃗₂)./r₂ + S⃗)

        if return_state_vector
		    return [r⃗₂, v⃗₂], get_orbital_elements(r⃗₂, v⃗₂, μ=μ, elements_type=elements_type)
        else
            return get_orbital_elements(r⃗₂, v⃗₂, μ=μ, elements_type=elements_type)
        end
	end
end