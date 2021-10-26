import numpy as np

import sys
sys.path.insert(1, '../Numerical-Methods/')
from solve_ode_ivp import rk_45


def periapsis_to_true_anomaly(eccentricity, angular_momentum, true_anomaly, mu=398600):
    orbit_period = (2 * np.pi / mu ** 2) * (angular_momentum / np.sqrt(1 - eccentricity ** 2)) ** 3

    eccentric_anomaly = 2 * np.arctan(np.sqrt((1 - eccentricity) / (1 + eccentricity) * np.tan(true_anomaly / 2)))
    mean_anomaly = eccentric_anomaly - eccentricity * np.sin(eccentric_anomaly)

    time = (mean_anomaly / (2 * np.pi)) * orbit_period

    return time


def two_body_propagator(t_init,
                        t_final,
                        mass_1,
                        mass_2,
                        initial_conditions,
                        grav_constant=6.67259e-11,
                        step_size=None,
                        tolerance=0.2,
                        beta=0.8):

    """
    Solves a general two-body problem using the rk_45 ode solver
    Initial conditions is a list of the form: [x_1, y_1, z_1, vx_1, vy_1, vz_1, x_2, y_2, z_2, vx_2, vy_2, vz_2]

    Returns a list, the first element being an array of each input variable, the second element being the times of each
        point. The solution array is the same form as the input
    """

    mu_1 = grav_constant * mass_1
    mu_2 = grav_constant * mass_2

    def differential_system(t, s):
        x_1, y_1, z_1, vx_1, vy_1, vz_1, x_2, y_2, z_2, vx_2, vy_2, vz_2 = s

        r = np.sqrt((x_2 - x_1) ** 2 + (y_2 - y_1) ** 2 + (z_2 - z_1) ** 2)

        return_vector = [
            vx_1,
            vy_1,
            vz_1,
            mu_2 * (x_2 - x_1) / r ** 3,
            mu_2 * (y_2 - y_1) / r ** 3,
            mu_2 * (z_2 - z_1) / r ** 3,
            vx_2,
            vy_2,
            vz_2,
            mu_1 * (x_1 - x_2) / r ** 3,
            mu_1 * (y_1 - y_2) / r ** 3,
            mu_1 * (z_1 - z_2) / r ** 3,
        ]

        return return_vector

    solution_out, time_out = rk_45(differential_system, initial_conditions, t_init, t_final, step_size, tolerance, beta)

    return solution_out, time_out


def three_body_cr_propagator(t_init,
                             t_final,
                             mass_1,
                             mass_2,
                             r_12,
                             initial_conditions,
                             grav_constant=6.67259e-11,
                             step_size=None,
                             tolerance=0.2,
                             beta=0.8):

    """
    Solves the circular restricted three-body problem using the rk_45 ode solver
    Initial conditions is a list of the form: [x, y, z, v_x, v_y, v_z]

    Returns a list, the first element being an array of each input variable, the second element being the times of each
        point. The solution array is the same form as the input
    """

    mu_1 = grav_constant * mass_1
    mu_2 = grav_constant * mass_2

    pi_1 = mass_1 / (mass_1 + mass_2)
    pi_2 = mass_2 / (mass_1 + mass_2)

    omega = np.sqrt(grav_constant * (mass_1 + mass_2) / r_12 ** 3)

    def differential_system(t, s):
        x, y, z, v_x, v_y, v_z = s

        r_1 = np.sqrt((x + pi_2 * r_12) ** 2 + y ** 2 + z ** 2)
        r_2 = np.sqrt((x - pi_1 * r_12) ** 2 + y ** 2 + z ** 2)

        return_vector = [
            v_x,
            v_y,
            v_z,
            2 * omega * v_y + x * omega ** 2 - (mu_1 / r_1 ** 3) * (x + pi_2 * r_12) - (mu_2 / r_2 ** 3) * (
                    x - pi_1 * r_12),
            -2 * omega * v_x + y * omega ** 2 - (mu_1 / r_1 ** 3) * y - (mu_2 / r_2 ** 3) * y,
            -(mu_1 / r_1 ** 3) * z - (mu_2 / r_2 ** 3) * z
        ]

        return return_vector

    solution_out, time_out = rk_45(differential_system, initial_conditions, t_init, t_final, step_size, tolerance, beta)

    return solution_out, time_out


def three_body_cr_thrust_propagator(t_init,
                                    t_final,
                                    initial_conditions,
                                    thrust_func,
                                    m_satellite,
                                    mass_1=5.974e24,
                                    mass_2=73.48e21,
                                    r_12=384400,
                                    grav_constant=6.6742e-11,
                                    return_thrust=False,
                                    step_size=None,
                                    tolerance=0.2,
                                    beta=0.8):

    """
    Solves the circular restricted three-body problem using the rk_45 ode solver, with included thrust terms to
    simulate spacecraft engines.

    Args:
        t_init: Starting time, in seconds
        t_final: Final time, in seconds
        initial_conditions: List of initial conditions of the spacecraft, in the form [x, y, z, v_x, v_y, v_z], in meters
        thrust_func: Function used to calculate the spacecraft's thrust, in N, defined as:
            thrust_func(t, x, y, z, v_x, v_y, v_z, m_satellite, omega, r_1, r_2, r_12, mu_1, mu_2, pi_1, pi_2)
        m_satellite: Mass of the spacecraft, in kg
        mass_1: Mass of the first body, in kg, mass of the Earth by default
        mass_2: Mass of the second body, in kg, mass of the Moon by default
        r_12: Distance between the two masses, distance between the Earth and the Moon by default
        grav_constant: The gravitational constant, 6.6742e-11 m^3 / k s^2 by default
        return_thrust: States whether to return the thrust components of the spacecraft along with the solution
        step_size: Initial step size to be input into the ODE solver
        tolerance: Tolerance of the ODE solver
        beta: Beta constant of the ODE solver

    Returns:
        solution_out: List of numpy arrays containing the variable solutions at each time step, in the form:
            [x, y, z, v_x, v_y, v_z]. A single array can be extracted with solution = solution_out.T[i]
        time_out: Numpy array of times
        thrust_out: List of numpy arrays containing the magnitude of thrust in each direction, not returned by default,
            the return_thrust argument must be set to True
    """

    mu_1 = grav_constant * mass_1
    mu_2 = grav_constant * mass_2

    pi_1 = mass_1 / (mass_1 + mass_2)
    pi_2 = mass_2 / (mass_1 + mass_2)

    omega = np.sqrt(grav_constant * (mass_1 + mass_2) / r_12 ** 3)

    def differential_system(t, s):
        x, y, z, v_x, v_y, v_z = s

        r_1 = np.sqrt((x + pi_2 * r_12) ** 2 + y ** 2 + z ** 2)
        r_2 = np.sqrt((x - pi_1 * r_12) ** 2 + y ** 2 + z ** 2)

        thrust_x, thrust_y, thrust_z = thrust_func(t, x, y, z, v_x, v_y, v_z, m_satellite, omega, r_1, r_2, r_12, mu_1, mu_2, pi_1, pi_2)

        return_vector = [
            v_x,
            v_y,
            v_z,
            thrust_x / m_satellite + 2 * omega * v_y + x * omega ** 2 - (mu_1 / r_1 ** 3) * (x + pi_2 * r_12) - (mu_2 / r_2 ** 3) * (x - pi_1 * r_12),
            thrust_y / m_satellite - 2 * omega * v_x + y * omega ** 2 - (mu_1 / r_1 ** 3) * y - (mu_2 / r_2 ** 3) * y,
            thrust_z / m_satellite - (mu_1 / r_1 ** 3) * z - (mu_2 / r_2 ** 3) * z
        ]

        return return_vector

    solution_out, time_out = rk_45(differential_system, initial_conditions, t_init, t_final, step_size, tolerance, beta)

    if return_thrust:
        # this is a hacky solution, other option is to write a custom rk_45 solver that appends each thrust value
        # when the solver approaches a solution at each time increment
        thrust_x_vec = []
        thrust_y_vec = []
        thrust_z_vec = []

        x_out = solution_out.T[0]
        y_out = solution_out.T[1]
        z_out = solution_out.T[2]

        vx_out = solution_out.T[3]
        vy_out = solution_out.T[4]
        vz_out = solution_out.T[5]

        for i in range(len(time_out)):
            t, x, y, z, v_x, v_y, v_z = time_out[i], x_out[i], y_out[i], 0, vx_out[i], vy_out[i], 0

            r_1 = np.sqrt((x + pi_2 * r_12) ** 2 + y ** 2 + z ** 2)
            r_2 = np.sqrt((x - pi_1 * r_12) ** 2 + y ** 2 + z ** 2)

            thrust_x, thrust_y, thrust_z = thrust_func(t, x, y, z, v_x, v_y, v_z, m_satellite, omega, r_1, r_2, r_12, mu_1, mu_2, pi_1, pi_2)

            thrust_x_vec.append(thrust_x)
            thrust_y_vec.append(thrust_y)
            thrust_z_vec.append(thrust_z)

        thrust_out = [thrust_x_vec, thrust_y_vec, thrust_z_vec]

        return solution_out, time_out, thrust_out

    else:
        return solution_out, time_out


def three_body_propagator(t_init,
                          t_final,
                          mass_1,
                          mass_2,
                          mass_3,
                          initial_conditions,
                          grav_constant=6.67259e-11,
                          step_size=None,
                          tolerance=0.2,
                          beta=0.8):

    """
    Solves the general three-body problem using the rk_45 ode solver
    Initial conditions is a list of the form:
        [x_1, y_1, z_1, vx_1, vy_1, vz_1, x_2, y_2, z_2, vx_2, vy_2, vz_2, x_3, y_3, z_3, vx_3, vy_3, vz_3]

    Returns a list, the first element being an array of each input variable, the second element being the times of each
        point. The solution array is the same form as the input
    """

    mu_1 = grav_constant * mass_1
    mu_2 = grav_constant * mass_2
    mu_3 = grav_constant * mass_3

    def differential_system(t, s):
        x_1, y_1, z_1, vx_1, vy_1, vz_1, x_2, y_2, z_2, vx_2, vy_2, vz_2, x_3, y_3, z_3, vx_3, vy_3, vz_3 = s

        r_12 = np.sqrt((x_2 - x_1) ** 2 + (y_2 - y_1) ** 2 + (z_2 - z_1) ** 2)
        r_13 = np.sqrt((x_3 - x_1) ** 2 + (y_3 - y_1) ** 2 + (z_3 - z_1) ** 2)
        r_32 = np.sqrt((x_3 - x_2) ** 2 + (y_3 - y_2) ** 2 + (z_3 - z_2) ** 2)

        return_vector = [
            vx_1,
            vy_1,
            vz_1,
            mu_2 * (x_2 - x_1) / r_12 ** 3 + mu_3 * (x_3 - x_1) / r_13 ** 3,
            mu_2 * (y_2 - y_1) / r_12 ** 3 + mu_3 * (y_3 - y_1) / r_13 ** 3,
            mu_2 * (z_2 - z_1) / r_12 ** 3 + mu_3 * (z_3 - z_1) / r_13 ** 3,
            vx_2,
            vy_2,
            vz_2,
            mu_1 * (x_1 - x_2) / r_12 ** 3 + mu_3 * (x_3 - x_2) / r_32 ** 3,
            mu_1 * (y_1 - y_2) / r_12 ** 3 + mu_3 * (y_3 - y_2) / r_32 ** 3,
            mu_1 * (z_1 - z_2) / r_12 ** 3 + mu_3 * (z_3 - z_2) / r_32 ** 3,
            vx_3,
            vy_3,
            vz_3,
            mu_1 * (x_1 - x_3) / r_13 ** 3 + mu_2 * (x_2 - x_3) / r_32 ** 3,
            mu_1 * (y_1 - y_3) / r_13 ** 3 + mu_2 * (y_2 - y_3) / r_32 ** 3,
            mu_1 * (z_1 - z_3) / r_13 ** 3 + mu_2 * (z_2 - z_3) / r_32 ** 3
        ]

        return return_vector

    solution_out, time_out = rk_45(differential_system, initial_conditions, t_init, t_final, step_size, tolerance, beta)

    return solution_out, time_out
