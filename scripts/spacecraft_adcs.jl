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