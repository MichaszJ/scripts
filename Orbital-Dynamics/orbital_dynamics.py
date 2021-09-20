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


def two_body_propagator(t_init, t_final, mass_1, mass_2, initial_conditions, grav_constant=6.67259e-11, step_size=None, tolerance=0.2, beta=0.8):
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
