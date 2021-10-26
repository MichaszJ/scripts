function two_body_propagator(t_init, t_final, mass_1, mass_2, initial_conditions; grav_constant=6.67259e-11, tolerance=0.2, beta=0.8)

    mu_1 = grav_constant * mass_1
    mu_2 = grav_constant * mass_2

    function differential_system(t, s)
        x_1, y_1, z_1, vx_1, vy_1, vz_1, x_2, y_2, z_2, vx_2, vy_2, vz_2 = s

        r = sqrt((x_2 - x_1)^2 + (y_2 - y_1)^2 + (z_2 - z_1)^2)

        return_vector = [
            vx_1,
            vy_1,
            vz_1,
            mu_2 * (x_2 - x_1) / r^3,
            mu_2 * (y_2 - y_1) / r^3,
            mu_2 * (z_2 - z_1) / r^3,
            vx_2,
            vy_2,
            vz_2,
            mu_1 * (x_1 - x_2) / r^3,
            mu_1 * (y_1 - y_2) / r^3,
            mu_1 * (z_1 - z_2) / r^3,
        ]

        return return_vector
    end

    solution_out, time_out = rk_45(differential_system, initial_conditions, t_init, t_final, tolerance, beta)

    return solution_out, time_out
end


function three_body_cr_propagator(t_init, t_final, mass_1, mass_2, r_12, initial_conditions; grav_constant=6.67259e-11, tolerance=0.2, beta=0.8)
    mu_1 = grav_constant * mass_1
    mu_2 = grav_constant * mass_2

    pi_1 = mass_1 / (mass_1 + mass_2)
    pi_2 = mass_2 / (mass_1 + mass_2)

    omega = sqrt(grav_constant * (mass_1 + mass_2) / r_12 ^ 3)

    function differential_system(t, s)
        x, y, z, v_x, v_y, v_z = s

        r_1 = sqrt((x + pi_2 * r_12) ^ 2 + y ^ 2 + z ^ 2)
        r_2 = sqrt((x - pi_1 * r_12) ^ 2 + y ^ 2 + z ^ 2)

        return_vector = [
            v_x,
            v_y,
            v_z,
            2 * omega * v_y + x * omega ^ 2 - (mu_1 / r_1 ^ 3) * (x + pi_2 * r_12) - (mu_2 / r_2 ^ 3) * (x - pi_1 * r_12),
            -2 * omega * v_x + y * omega ^ 2 - (mu_1 / r_1 ^ 3) * y - (mu_2 / r_2 ^ 3) * y,
            -(mu_1 / r_1 ^ 3) * z - (mu_2 / r_2 ^ 3) * z
        ]

        return return_vector
    end

    solution_out, time_out = rk_45(differential_system, initial_conditions, t_init, t_final, tolerance, beta)

    return solution_out, time_out
end