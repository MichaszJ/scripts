include("solve_ode.jl")

function rocket_launch_trajectory(rocket_params, initial_conditions, t0, tf, thrust_func, ρ_func, D_func, g_func, m_func, h_pitchover; r_planet=6378.0e3, tolerance=50.0, beta=0.8, solver_method="rk_45")
    # defining the system of differential equations
    function differential_system(t::Float64, s0::Vector{Float64})
        # extracting values from input
        v, γ, h, x = s0

        # getting physical/structural rocket parameters
        # TODO: expand and add more parameters
        launch_mass, mass_ratio, drag_coeff, frontal_area = rocket_params
        final_mass = launch_mass / mass_ratio

        # calculating values from functions
        density = ρ_func(h)
        D = D_func(density, v)
        g = g_func(h)
        m = m_func(t)

        

        # getting thrust from universal thrust function
        m > final_mass ? thrust = thrust_func(t, v, γ, h, x, density, D, g, m) : thrust = 0

        # TODO: implement returning stage in rocket flight, e.g., vertical flight, gravity turn, MECO, etc.

        # returning differentials based on whether rocket is in initial burn or gravity turn
        # TODO: implement functioanlity for multiple stages
        if h <= h_pitchover
            return [
                thrust/m - D/m - g,
                0,
                v,
                0,
            ]
        elseif h > h_pitchover && m > final_mass
            return [
                thrust/m - D/m - g*sin(γ),
                -(1/v) * (g - v^2 / (r_planet + h)) * cos(γ),
                v * sin(γ),
                (r_planet / (r_planet + h)) * v * cos(γ),
            ]
        else
            return [
                -D/m - g*sin(γ),
                -(1/v) * (g - v^2 / (r_planet + h)) * cos(γ),
                v * sin(γ),
                (r_planet / (r_planet + h)) * v * cos(γ),
            ]
        end
    end

    s_out, t_out = rk_solver(differential_system, initial_conditions, t0, tf, solver_method=solver_method, tolerance=tolerance, beta=beta)
    
    return s_out, t_out
end

function diffeq_nonimpulsive_orbital_maneuver(initial_conditions, time_span, params; solver_args...)
	function differential_system!(du, u, p, t)
		# unpacking initial conditions and parameters
		x, y, z, ẋ, ẏ, ż, m = u
		μ, Isp, Tx_func, Ty_func, Tz_func = p
	
		# defining useful values
		Tx = Tx_func(t, x, y, z, ẋ, ẏ, ż)
		Ty = Ty_func(t, x, y, z, ẋ, ẏ, ż)
		Tz = Tz_func(t, x, y, z, ẋ, ẏ, ż)
		r = sqrt(x^2 + y^2 + z^2)
		v = sqrt(ẋ^2 + ẏ^2 + ż^2)
		
		# defining system of DEs
		du[1] = ẋ
		du[2] = ẏ
		du[3] = ż
		du[4] = -μ*x / r^3 + (Tx / m) * (ẋ / v)
		du[5] = -μ*y / r^3 + (Ty / m) * (ẏ / v)
		du[6] = -μ*z / r^3 + (Tz / m) * (ż / v)
		du[7] = -1e3 * sqrt(Tx^2 + Ty^2 + Tz^2) / (Isp * 9.80665)
	end

	problem = ODEProblem(differential_system!, initial_conditions, time_span, params)
	solution = solve(problem; solver_args...)
	
	return solution
end