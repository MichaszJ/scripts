include("solve_ode_ivp.jl")

function rocket_launch_trajectory(rocket_params, initial_conditions, t0, tf, ρ_func, D_func, g_func, m_func, h_pitchover; r_planet=6378.0e3, tolerance=50.0, beta=0.8, solver_method="rk_45")
    # defining the system of differential equations
    function differential_system(t::Float64, s0::Vector{Float64})
        # extracting values from input
        v, γ, h, x = s0

        # getting rocket parameters
        # TODO: expand and add more parameters
        launch_mass, mass_ratio, thrust, mass_flux, drag_coeff, frontal_area = rocket_params
        
        # calculating values from functions
        density = ρ_func(h)
        D = D_func(density, v)
        g = g_func(h)
        m = m_func(t)

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
        elseif h > h_pitchover && m > launch_mass/mass_ratio
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