function par_1(θ)
    return [
        1 0 0
        0 cos(θ) sin(θ)
        0 -sin(θ) cos(θ)
    ]
end

function par_2(θ)
     return [
        cos(θ) 0 -sin(θ)
        0 1 0
        sin(θ) 0 cos(θ)
    ]
end

function par_3(θ)
    return [
        cos(θ) sin(θ) 0 
        -sin(θ) cos(θ) 0 
        0 0 1
    ]
end

function cross_mat(v)
    @assert length(v) == 3 "Vectors of length 3 required"
    return [
        0 -v[3] v[2]
        v[3] 0 -v[1]
        -v[2] v[1] 0
    ]
end

function trace(matrix)
	@assert size(matrix,1) == size(matrix,2) "Input must be a square matrix" 
	return sum(matrix[i,i] for i in 1:size(matrix,1))
end

# convert quaternion to direction cosine matrix
function q2DCM(q)
    q1, q2, q3, q4 = q
    return [
        q4^2+q1^2-q2^2-q3^2 2*(q1*q2 + q4*q3) 2*(q1*q3 - q4*q2)
        2*(q1*q2 - q4*q3) (q4^2-q1^2+q2^2-q3^2) 2*(q2*q3+q4*q1)
        2*(q1*q3 + q4*q2) 2*(q2*q3-q4*q1) (q4^2-q1^2-q2^2+q3^2)
    ]
end

# convert classic rodrigues parameter to direction cosine matrix
function p2DCM(p)
    q4 = 1 / sqrt(1 + norm(p)^2)
    qv = p ./ sqrt(1 + norm(p)^2)
    
    return q2DCM([qv[1], qv[2], qv[3], q4])
end

function ecef2geocentric_latlon(x, y, z)
    λ = atan(y, x)
    p = sqrt(x^2 + y^2)
    ϕc = atan(p, z)

    return π/2 - ϕc, λ 
end

function ecef2geodetic_latlon(x, y, z; iters=5)
    λ = atan(y, x)
    
    r = sqrt(x^2 + y^2 + z^2)
    p = sqrt(x^2 + y^2)
    
    ϕc = atan(p, z)
    ϕ = ϕc
    
    a = 6378.137
    esq = 6.69437999014e-3
    for _ in 1:iters
        Rn = a / sqrt(1 - esq * sin(ϕ)^2)
        h = p/cos(ϕ) - Rn
        ϕ = atan((z/p) / (1 - esq * (Rn / (Rn + h))))
    end
    
    return π/2 - ϕc, λ 
end

function local_polar(R)
    x, y, z = R
    δ, λ = ecef2geocentric_latlon(x, y, z)
    
    R_BE = par_3((2*π - λ)) * par_2(δ-π/2)
    
    x̂B = R_BE * [1, 0, 0]
    ŷB = R_BE * [0, 1, 0]
    ẑB = R_BE * [0, 0, 1]
    
    return x̂B, ŷB, ẑB
end

function polar2ecef2(R, v1, v2, v3)
    x, y, z = R
    δ, λ = ecef2geocentric_latlon(x, y, z)
    
    R_EB = par_2(π/2 - δ) * par_3(λ - 2*π)
    
    ex = [(R_EB * v1)[1], 0.0, 0.0]
    ey = [0.0, (R_EB * v2)[2], 0.0]
    ez = [0.0, 0.0, (R_EB * v3)[3]]
    
    return ex, ey, ez
end

function C_LVLH_ECEF(Rx, Ry, Rz, Vx, Vy, Vz)
    o3 = [Rx, Ry, Rz] ./ norm([Rx, Ry, Rz])
    o2 = cross([Rx, Ry, Rz], [Vx, Vy, Vz]) ./ norm(cross([Rx, Ry, Rz], [Vx, Vy, Vz]))
    o1 = cross(o2, o3)

    return [o1 o2 o3]
end

function dcm_to_axisangle(DCM)
    ϕ = 0.5 * (trace(DCM) - 1)
    a_vec = [
        DCM[2,3] - DCM[3,2]
        DCM[3,1] - DCM[1,3]
        DCM[1,2] - DCM[2,1]
    ] ./ (2 * sin(ϕ))
    
    return ϕ, a_vec
end

function dcm_to_eulerparameter(DCM)
    ϕ, a_vec = dcm_to_axisangle(DCM)
    
    η = cos(ϕ / 2)
    ϵ = sin(ϕ / 2) .* a_vec
    
    return η, ϵ
end

function dcm_to_quaternion(DCM)
    η, ϵ = dcm_to_eulerparameter(DCM)
    
    return [
        η
        ϵ
    ]
end

function quest_method(ref_vecs, mea_vecs; iter=5, weights=nothing)
	if weights === nothing
		weights = ones(length(ref_vecs))
	end

	B = sum(weights .* mea_vecs .* transpose.(ref_vecs))

	S = B + B'

	σ = trace(B)

	Z = [
        B[2,3] - B[3,2]
        B[3,1] - B[1,3]
        B[1,2] - B[2,1]
    ]

	K = [
        S-σ*I(3) Z
        Z' σ
    ]

	λ = [sum(weights) for i in 1:4]

	a = σ^2 - trace(adjoint(S))
	b = σ^2 + Z' * Z
	c = 8*det(B)
	d = Z' * S * S * Z

	_poly(λ) = (λ.^2 .- a) .* (λ.^2 .- b) .- c.*λ .+ (c*σ - d)
	_poly_prime(λ) = 4 .* λ .^3 .- 2*(a + b) .* λ .- c

	for _ in 1:iter
        λ = λ - _poly(λ) ./ _poly_prime(λ)
    end

	λmax = maximum(λ)
    p = inv((λmax + σ)*I(3) - S) * Z
    
    q̄ = (1 / sqrt(1 + p' * p)) .* [p; 1]
    
    return q̄
end

function q_method(ref_vecs, mea_vecs; weights=nothing)
	if weights === nothing
		weights = ones(length(ref_vecs))
	end

	B = sum(weights .* mea_vecs .* transpose.(ref_vecs))

	S = B + B'

	σ = trace(B)

	Z = [
        B[2,3] - B[3,2]
        B[3,1] - B[1,3]
        B[1,2] - B[2,1]
    ]

	K = [
        S-σ*I(3) Z
        Z' σ
    ]

	λmax = maximum(eigvals(K))
	
    p = inv((λmax + σ)*I(3) - S) * Z
    
    q̄ = (1 / sqrt(1 + p' * p)) .* [p; 1]
    
    return q̄
end

function diffeq_euler_quat(initial_conditions, time_span, params; solver_args...)
    function _differential_system(u, p, t)
        q0, q1, q2, q3, ω1, ω2, ω3 = u
        J1, J2, J3, M1, M2, M3 = p

        return SA[
            -0.5 * (q1*ω1 + q2*ω2 + q3*ω3),
            0.5 * (q0*ω1 + q2*ω3 - q3*ω2),
            0.5 * (q0*ω2 - q1*ω3 + q3*ω1),
            0.5 * (q0*ω3 + q1*ω2 - q2*ω1),
            (M1(u, p, t) + (J2 - J3) * ω2 * ω3) / J1,
            (M2(u, p, t) + (J3 - J1) * ω1 * ω3) / J2,
            (M3(u, p, t) + (J1 - J2) * ω1 * ω2) / J3
        ]
    end

    problem = ODEProblem(_differential_system, initial_conditions, time_span, params)
    solution = solve(problem; solver_args...)

    return solution
end

function diffeq_euler_eulerangle(initial_conditions, time_span, params; solver_args...)
    function _differential_system(u, p, t)
        ϕ, θ, ψ, ω_x, ω_y, ω_z = u
        J_x, J_y, J_z, M_x, M_y, M_z = p
        
        return SA[
            ω_x + ω_z * tan(θ)*cos(ϕ) + ω_y*tan(θ)*sin(ϕ),
            ω_y*cos(ϕ) - ω_z*sin(ϕ),
            ω_z*sec(θ)*cos(ϕ) + ω_y*sec(θ)*sin(ϕ),
            (M_x(u, p, t) + (J_y - J_z) * ω_y * ω_z) / J_x,
            (M_y(u, p, t) + (J_z - J_x) * ω_x * ω_z) / J_y,
            (M_z(u, p, t) + (J_x - J_y) * ω_x * ω_y) / J_z
        ]
    end

    problem = ODEProblem(_differential_system, initial_conditions, time_span, params)
    solution = solve(problem; solver_args...)

    return solution
end