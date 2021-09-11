import numpy as np
import matplotlib.pyplot as plt

from numerical_methods import find_root


def y_func(xin):
    return (xin - 2) ** 2


x = np.linspace(0, 5, 100)
y = y_func(x)

root, xs, ys, errors = find_root(x, y, 1, epsilon=1e-5, max_iterations=250)
print(root)

plt.semilogy(errors)
plt.show()

