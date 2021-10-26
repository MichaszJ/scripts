import numpy as np


# basic 4th order runge-kutta method
def rk_4(ode_func, y_init, t_init, t_final, step_size=None):
    if step_size is None:
        step_size = (t_final - t_init) / 100

    t_out = [t_init]
    y_out = [np.array(y_init)]

    while t_init < t_final:
        f1 = np.array(ode_func(t_init, y_init))
        f2 = np.array(ode_func(t_init + 0.5 * step_size, y_init + 0.5 * step_size * f1))
        f3 = np.array(ode_func(t_init + 0.5 * step_size, y_init + 0.5 * step_size * f2))
        f4 = np.array(ode_func(t_init + step_size, y_init + step_size * f3))

        y_init = y_init + step_size * ((1 / 6) * f1 + (1 / 3) * f2 + (1 / 3) * f3 + (1 / 6) * f4)
        t_init = t_init + step_size

        y_out.append(y_init)
        t_out.append(t_init)

    return [np.array(y_out), np.array(t_out)]


# runge-kutta with variable step size
def rk_45(ode_func, y_init, t_init, t_final, step_size=None, tolerance=0.2, beta=0.8):
    if step_size is None:
        step_size = (t_final - t_init) / 100

    a = [0, 1 / 4, 3 / 8, 12 / 13, 1, 1 / 2]
    b = [
        [0, 0, 0, 0, 0],
        [1 / 4, 0, 0, 0, 0],
        [3 / 32, 9 / 32, 0, 0, 0],
        [1932 / 2197, -7200 / 2197, 7296 / 2197, 0, 0],
        [439 / 216, -8, 3680 / 513, -845 / 4104, 0],
        [-8 / 27, 2, -3544 / 2565, 1859 / 4104, -11 / 40]
    ]
    cs = [25 / 216, 0, 1408 / 2565, 2197 / 4104, -1 / 5, 0]
    c = [16 / 135, 0, 6656 / 12825, 28561 / 56430, -9 / 50, 2 / 55]

    t_out = [t_init]
    y_out = [np.array(y_init)]

    f_vector = [0] * 6

    while t_init < t_final:
        for i in range(len(f_vector)):
            t_delta = a[i] * step_size
            y_delta = step_size * np.array(sum([b[i][j] * f_vector[j] for j in range(i)]))

            f_vector[i] = np.array(ode_func(t_init + t_delta, y_init + y_delta))

        truncation_vector = np.array([abs(step_size * (c[i] - cs[i]) * f_vector[i]) for i in range(len(f_vector))])
        truncation_error = truncation_vector.max()

        if truncation_error > tolerance:
            step_size = step_size * beta * np.power((tolerance / truncation_error), 1 / 5)
        else:
            y_init = y_init + step_size * sum([c[i] * f_vector[i] for i in range(len(c))])
            t_init = t_init + step_size

            y_out.append(y_init)
            t_out.append(t_init)

            step_size = step_size * beta * np.power((tolerance / truncation_error), 1 / 5)

    return [np.array(y_out), np.array(t_out)]
