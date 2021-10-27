include("solve_ode_ivp.jl")

function rocket_launch_trajectory(rocket_params, initial_conditions, t0, tf, ρ_func, D_func, g_func, m_func, h_pitchover; r_planet=6378.0e3, tolerance=50.0, beta=0.8, solver_method="rk_45")
    # defining the system of differential equations
    function differential_system(t::Float64, s0::Vector{Float64})
        # extracting values from input
        v, γ, h, x = s0

        launch_mass, mass_ratio, thrust, mass_flux, drag_coeff, frontal_area = rocket_params
        
        # calculating values from functions
        density = density_func(h)
        D = drag_func(density, v)
        g = gravity_func(h)
        m = mass_func(t)

        # returning differentials based on whether rocket is in initial burn or gravity turn
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

    if solver_method == "rk_45"
        s_out, t_out = rk_45(differential_system, initial_conditions, t0, tf, tolerance=tolerance, beta=beta)
    elseif solver_method == "rk_56"
        s_out, t_out = rk_56(differential_system, initial_conditions, t0, tf, tolerance=tolerance, beta=beta)
    elseif solver_method == "rk_78"
        s_out, t_out = rk_78(differential_system, initial_conditions, t0, tf, tolerance=tolerance, beta=beta)
    else
        println("Invalid solver method")
        s_out, t_out = nothing, nothing
    end
    return s_out, t_out
end