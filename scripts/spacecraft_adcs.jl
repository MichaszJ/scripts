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