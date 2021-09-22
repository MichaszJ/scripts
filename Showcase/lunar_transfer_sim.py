import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

from matplotlib import animation
from matplotlib.animation import PillowWriter

import sys

sys.path.insert(1, '../Numerical-Methods/')
from solve_ode_ivp import rk_45

sys.path.insert(1, '../Utilities/')
from animations import constant_time_transform

sns.set()

# defining custom propagator
def three_body_cr_thrust_propagator(t_init,
                                    t_final,
                                    mass_1,
                                    mass_2,
                                    r_12,
                                    initial_conditions,
                                    thrust_func,
                                    m_satellite,
                                    grav_constant=6.67259e-11,
                                    step_size=None,
                                    tolerance=0.2,
                                    beta=0.8):
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
            thrust_y / m_satellite + -2 * omega * v_x + y * omega ** 2 - (mu_1 / r_1 ** 3) * y - (mu_2 / r_2 ** 3) * y,
            thrust_z / m_satellite + -(mu_1 / r_1 ** 3) * z - (mu_2 / r_2 ** 3) * z
        ]

        return return_vector

    solution_out, time_out = rk_45(differential_system, initial_conditions, t_init, t_final, step_size, tolerance, beta)

    return solution_out, time_out


# defining thrust function
def thrust(t, x, y, z, v_x, v_y, v_z, m_satellite, omega, r_1, r_2, r_12, mu_1, mu_2, pi_1, pi_2):
    if (-100 <= y <= 100) and (1 <= np.sqrt(v_x ** 2 + v_y ** 2 + v_z ** 2) <= 2) and t < 300000:
        thrust_y = 5
    else:
        thrust_y = 0

    thrust_x = 0
    thrust_z = 0

    return thrust_x, thrust_y, thrust_z


# defining initial conditions, units in km, km/s
x = -4671
y = -6378 - 200
z = 0

v_x = 10.9148 * np.cos(np.radians(19))
v_y = -10.9148 * np.sin(np.radians(19))
v_z = 0

s0 = [x, y, z, v_x, v_y, v_z]


# running propagator
s_out_3b, t_out_3b = three_body_cr_thrust_propagator(t_init=0,
                                                     t_final=8 * 24 * 60 * 60,
                                                     mass_1=5.974e24,
                                                     mass_2=73.48e21,
                                                     r_12=384400,
                                                     initial_conditions=s0,
                                                     thrust_func=thrust,
                                                     m_satellite=1000,
                                                     grav_constant=6.6759e-20,
                                                     tolerance=5)

# extracting desired results
x_out = s_out_3b.T[0]
y_out = s_out_3b.T[1]


# converting results to constant time intervals, for use in animation
t_new = np.linspace(0, 8 * 24 * 60 * 60, 1000)
x_points, y_points = constant_time_transform([x_out, y_out], t_out_3b, t_new)


# creating animation, same animation featured in scripts/README
fig, ax = plt.subplots(figsize=(10, 5), dpi=150)

time_text = ax.text(0.87, 0.95, '', transform=ax.transAxes)
point = ax.scatter(x_points[0], y_points[0], s=2)


def animate(i):
    ax.plot(x_points[0:i], y_points[0:i], linewidth=0.75, c='b')
    point.set_offsets([x_points[i], y_points[i]])
    time_text.set_text(f't = {round(t_new[i] / (24 * 60 * 60), 2)} Days')


ax.set_title('Earth to Moon Transfer')

earth = plt.Circle((-4671, 0), 6378)
ax.add_artist(earth)

moon = plt.Circle((-4671 + 384400, 0), 1737)
ax.add_artist(moon)

ax.set_xlabel('x (km)')
ax.set_ylabel('y (km)')

ax.set_xlim(-15000, 400000)
ax.set_ylim(-20000, 100000)

plt.gca().set_aspect('equal')

ani = animation.FuncAnimation(fig, animate, frames=500)
ani.save('lunar_transfer.gif', writer='pillow', fps=30)
