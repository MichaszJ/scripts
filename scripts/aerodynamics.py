import numpy as np
import matplotlib.pyplot as plt
from skaero.atmosphere import coesa


def flight_envelope_graph(n_lim, wing_area, max_takeoff, cruise_mach, cruise_alt, max_mach, max_cl):
    alt_temperature = coesa.table(cruise_alt)[1]
    alt_density = coesa.table(cruise_alt)[3]
    sl_density = coesa.table(0)[3]

    eq_cruise_vel = np.sqrt(alt_density / sl_density) * cruise_mach * np.sqrt(1.4 * 287 * alt_temperature)
    eq_max_vel = np.sqrt(alt_density / sl_density) * max_mach * np.sqrt(1.4 * 287 * alt_temperature)

    n_vec_pos = np.linspace(0, n_lim)
    n_vec_neg = np.linspace(0, -1)

    def velocity(n):
        return np.sqrt((2 * abs(n) * (max_takeoff * 9.81)) / (sl_density * wing_area * max_cl))

    vel_pos = velocity(n_vec_pos)
    vel_neg = velocity(n_vec_neg)

    # plotting stall curves
    plt.plot(vel_pos, n_vec_pos, color='b', linewidth=0.75)
    plt.plot(vel_neg, n_vec_neg, color='b', linewidth=0.75)

    # plotting and annotating definite points
    plt.scatter(eq_cruise_vel, n_lim, color='b')
    plt.annotate(
        str(np.round(eq_cruise_vel, decimals=4)),
        xy=(eq_cruise_vel, n_lim), xycoords='data',
        xytext=(eq_cruise_vel - 6.5, n_lim + 0.4), textcoords='data',
        arrowprops=dict(arrowstyle="->", connectionstyle="arc3"),
    )

    plt.scatter(eq_max_vel, n_lim, color='b')
    plt.annotate(
        str(np.round(eq_max_vel, decimals=4)),
        xy=(eq_max_vel, n_lim), xycoords='data',
        xytext=(eq_max_vel + 10, n_lim), textcoords='data',
        arrowprops=dict(arrowstyle="->", connectionstyle="arc3"),
    )

    plt.scatter(vel_pos[-1], n_lim, color='b')
    plt.annotate(
        str(np.round(vel_pos[-1], decimals=4)),
        xy=(vel_pos[-1], n_lim), xycoords='data',
        xytext=(vel_pos[-1] - 25, n_lim), textcoords='data',
        arrowprops=dict(arrowstyle="->", connectionstyle="arc3"),
    )

    plt.scatter(vel_neg[-1], -1, color='b')
    plt.annotate(
        str(np.round(vel_neg[-1], decimals=4)),
        xy=(vel_neg[-1], -1), xycoords='data',
        xytext=(vel_neg[-1], -0.5), textcoords='data',
        arrowprops=dict(arrowstyle="->", connectionstyle="arc3"),
    )

    plt.scatter(eq_cruise_vel, -1, color='b')

    # connecting points together
    plt.plot([vel_pos[-1], eq_max_vel, eq_max_vel, eq_cruise_vel, vel_neg[-1]], [n_lim, n_lim, 0, -1, -1], color='b',
             linewidth=0.75)
    plt.plot([eq_cruise_vel, eq_cruise_vel], [n_lim, -1], linestyle='--', color='orange', linewidth=0.75)

    plt.xlabel('V (EAS) [m/s]')
    plt.ylabel('n')
