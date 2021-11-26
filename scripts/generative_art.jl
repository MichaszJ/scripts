function create_spirograph(R, r, ρ, t)
	k = r/R
	l = ρ/r

	x = @. R * ((1 - k) * cos(t) + l*k*cos(t * (1 - k)/ k))
	y = @. R * ((1 - k) * sin(t) - l*k*sin(t * (1 - k)/ k))

    return x, y
end