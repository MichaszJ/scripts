import numpy as np

from numerical_methods import find_root_data, find_root_function


def y_func(xin):
    return xin - 3


x = np.linspace(0, 5, 500)
y = y_func(x)

root1 = find_root_data(x, y, [1, 4.5], method='bisection')
print(root1)

root2 = find_root_data(x, y, 4, method='secant')
print(root2)

root3 = find_root_function(y_func, [1, 4.5], method='bisection')
print(root3)

root4 = find_root_function(y_func, 2, method='secant')
print(root4)

