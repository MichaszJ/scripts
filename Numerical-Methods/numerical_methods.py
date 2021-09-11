import numpy as np


def find_root(x, y, initial_guess, epsilon=1e-4, method='newton', max_iterations=100, debug_text=False):
    x_final = initial_guess
    errors = []
    xs = []
    ys = []

    def interpolate(x_in, x_0, x_1, y_0, y_1):
        return y_0 + (x_in - x_0) * ((y_1 - y_0) / (x_1 - x_0))

    def error(value_current, value_previous):
        if value_previous is not None:
            err = abs((value_previous - value_current) / value_current)
        else:
            err = None

        return err

    if method == 'newton':
        diff_array = np.diff(x) / np.diff(y)

        iteration = 0
        current_error = 100

        x_n = initial_guess
        x_prev = None
        x_final = initial_guess

        while current_error > epsilon and iteration < max_iterations:
            # find closest value to initial guess in the given data
            for i, xi in enumerate(x):
                if xi > x_n > x[i - 1] or x_n == x[i]:
                    x1 = xi
                    x0 = x[i - 1]

                    y1 = y[i]
                    y0 = y[i - 1]

                    yp1 = diff_array[i]
                    yp0 = diff_array[i - 1]

                    break

            diff = interpolate(x_n, x0, x1, yp0, yp1)
            y_n = interpolate(x_n, x0, x1, y0, y1)
            iteration_error = error(x_n, x_prev)

            errors.append(iteration_error)
            xs.append(x_n)
            ys.append(y_n)

            if debug_text is True:
                print(f'Iteration: {iteration + 1} | Current x: {x_n} | Current y: {y_n} | Error: {iteration_error}')

            if iteration_error is not None and epsilon > iteration_error > -epsilon:
                break

            x_prev = x_n
            x_n = x_n - y_n / diff
            x_final = x_n

            iteration += 1

    return x_final, xs, ys, errors
