import numpy as np
import scipy.integrate as integrate
import scipy.interpolate as interpolate

def corrected_cl(cp_upper, cp_lower, mach_num, num_points=500):
    x_vec = np.linspace(0, 1, num_points)

    upper = cp_upper(x_vec)
    lower = cp_lower(x_vec)

    pg_correction = 1 / np.sqrt(1 - np.power(mach_num, 2))

    cp_difference = interpolate.UnivariateSpline(x_vec, upper - lower)
    lift_coefficient = -integrate.quad(cp_difference, 0, 1)[0]
    corrected_cl = lift_coefficient * pg_correction

    return corrected_cl

def corrected_distribution(cp_upper, cp_lower, mach_num, num_points=500):
    x_vec = np.linspace(0, 1, num_points)

    upper = cp_upper(x_vec)
    lower = cp_lower(x_vec)

    pg_correction = 1 / np.sqrt(1 - np.power(mach_num, 2))

    upper_corrected = upper * pg_correction
    lower_corrected = lower * pg_correction

    return [upper_corrected, lower_corrected]

def critical_pressure_coef(mach_num):
    k1 = 2 / (1.4 * np.power(mach_num, 2))
    k2 = (1 + 0.5 * (1.4 - 1) * np.power(mach_num, 2)) / (1 + 0.5 * (1.4 - 1))
    k3 = 1.4 / (1.4 - 1)

    cp_crit = k1 * ((np.power(k2, k3)) - 1)

    return cp_crit