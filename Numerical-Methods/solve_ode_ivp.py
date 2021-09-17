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
    ys = np.array(y_init)

    while t_init < t_final:
        f1 = np.array(ode_func(t_init, y_init))
        f2 = np.array(ode_func(t_init + a[1] * step_size, y_init + step_size * b[1][0] * f1))
        f3 = np.array(ode_func(t_init + a[2] * step_size, y_init + step_size * (b[2][0] * f1 + b[2][1] * f2)))
        f4 = np.array(
            ode_func(t_init + a[3] * step_size, y_init + step_size * (b[3][0] * f1 + b[3][1] * f2 + b[3][2] * f3)))
        f5 = np.array(ode_func(t_init + a[4] * step_size,
                               y_init + step_size * (b[4][0] * f1 + b[4][1] * f2 + b[4][2] * f3 + b[4][3] * f4)))
        f6 = np.array(ode_func(t_init + a[5] * step_size, y_init + step_size * (
                    b[5][0] * f1 + b[5][1] * f2 + b[5][2] * f3 + b[5][3] * f4 + b[5][4] * f5)))

        f_vector = [f1, f2, f3, f4, f5, f6]

        truncation_vector = np.array([abs(step_size * (c[i] - cs[i]) * f_vector[i]) for i in range(len(f_vector))])
        truncation_error = truncation_vector.max()

        if truncation_error > tolerance:
            step_size = step_size * beta * np.power((tolerance / truncation_error), 1 / 5)
        else:
            y_init = y_init + step_size * (c[0] * f1 + c[1] * f2 + c[2] * f3 + c[3] * f4 + c[4] * f5 + c[5] * f6)
            t_init = t_init + step_size

            y_out.append(y_init)
            t_out.append(t_init)

            step_size = step_size * beta * np.power((tolerance / truncation_error), 1 / 5)

    return [np.array(y_out), np.array(t_out)]
