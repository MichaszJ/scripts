import numpy as np


def find_root_data(x, y, initial_guess, epsilon=1e-4, method='secant', max_iterations=50, debug_text=False):
    x_final = initial_guess

    def interpolate(x_in, x_0, x_1, y_0, y_1):
        return y_0 + (x_in - x_0) * ((y_1 - y_0) / (x_1 - x_0))

    def relative_error(value_current, value_previous):
        if value_previous is not None:
            err = abs((value_previous - value_current) / value_current)
        else:
            err = None

        return err

    if method == 'secant':
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
            iteration_error = relative_error(x_n, x_prev)

            if debug_text is True:
                print(f'Iteration: {iteration + 1} | Current x: {x_n} | Current y: {y_n} | Error: {iteration_error}')

            if iteration_error is not None and epsilon > iteration_error > -epsilon:
                break

            x_prev = x_n
            x_n = x_n - y_n / diff
            x_final = x_n

            iteration += 1

    elif method == 'bisection':
        endpoints = initial_guess
        iteration = 0
        current_error = 100
        x_prev = None
        x_final = None

        while current_error > epsilon and iteration < max_iterations:
            x_a = endpoints[0]
            x_b = endpoints[1]
            x_midpoint = (x_a + x_b) / 2

            for i, xi in enumerate(x):
                if xi > x_midpoint > x[i - 1] or x_midpoint == x[i]:
                    y_midpoint = interpolate(x_midpoint, x[i - 1], x[i], y[i - 1], y[i])

            x_final = x_midpoint

            if debug_text is True:
                print(f'Iteration: {iteration + 1} | Current x: {x_midpoint} | Current y: {y_midpoint}')

            if abs(y_midpoint) < epsilon:
                break

            if y_midpoint < 0:
                endpoints = [x_midpoint, x_b]
            else:
                endpoints = [x_a, x_midpoint]

            iteration += 1

    else:
        print('Please enter a valid method')

    return x_final


def find_root_function(function, initial_guess, differential=1e-5, epsilon=1e-4, method='secant', max_iterations=50,
                       debug_text=False):
    x_final = initial_guess

    def relative_error(value_current, value_previous):
        if value_previous is not None:
            err = abs((value_previous - value_current) / value_current)
        else:
            err = None

        return err

    if method == 'secant':
        iteration = 0
        current_error = 100

        x_n = initial_guess
        x_prev = None
        x_final = initial_guess

        while current_error > epsilon and iteration < max_iterations:
            diff = (function(x_n + differential) - function(x_n - differential)) / differential
            iteration_error = relative_error(x_n, x_prev)

            if debug_text is True:
                print(f'Iteration: {iteration + 1} | Current x: {x_n} | Current y: {y_n} | Error: {iteration_error}')

            if iteration_error is not None and epsilon > iteration_error > -epsilon:
                break

            x_prev = x_n
            x_n = x_n - function(x_n) / diff
            x_final = x_n

            iteration += 1

    elif method == 'bisection':
        endpoints = initial_guess
        iteration = 0
        current_error = 100

        while current_error > epsilon and iteration < max_iterations:
            x_a = endpoints[0]
            x_b = endpoints[1]
            x_midpoint = (x_a + x_b) / 2

            y_midpoint = function(x_midpoint)

            x_final = x_midpoint

            if debug_text is True:
                print(f'Iteration: {iteration + 1} | Current x: {x_midpoint} | Current y: {y_midpoint}')

            if abs(y_midpoint) < epsilon:
                break

            if y_midpoint < 0:
                endpoints = [x_midpoint, x_b]
            else:
                endpoints = [x_a, x_midpoint]

            iteration += 1

    else:
        print('Please enter a valid method')

    return x_final
