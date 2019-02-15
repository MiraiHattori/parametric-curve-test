import numpy as np
import matplotlib.pyplot as plt

def calc_bspline(x, j, N, M, x_max, x_min):
    return calc_bspline_coeff_vector_base(x, j, N, N, M, x_max, x_min)

def calc_bspline_coeff_vector_base(x, j, N, _n, M, x_max, x_min):
    h = (x_max - x_min) * 1.0 / (M - N)
    xj = h * j + (M * x_min - N * x_max) * 1.0 / (M - N)
    xj_1 = xj + h
    xj_n_1 = xj_1 + _n * h
    if x < xj or xj_n_1 <= x:
        return np.zeros(_n + 1)
    elif _n == 0:
        return np.ones(1)
    else:
        v1 = (x - xj) * calc_bspline_coeff_vector_base(x, j, N, _n - 1, M, x_max, x_min)
        v2 = (xj_n_1 - x) * calc_bspline_coeff_vector_base(x, j + 1, N, _n - 1, M, x_max, x_min)
        return (np.hstack((v1 - v2, np.zeros(1))) + np.hstack((np.zeros(1), xj_n_1 * v2 - xj * v1))) * 1.0 / (_n * h)


N = 2
M = 5
x_max = 5
x_min = 0

h = (x_max - x_min) * 1.0 / (M - N)

print(calc_bspline(0, 0, N, M, x_max, x_min))

"""
fig, ax = plt.subplots()
xx = np.linspace(-1, 8, 100)
c = [1.0, 1.0, 1.0]
ax.plot(xx, [calc_bspline(x, 0, N, M, x_max, x_min) * c[0] for x in xx], '-', lw=2, label='bspline-base-0')
ax.plot(xx, [calc_bspline(x, 1, N, M, x_max, x_min) * c[1] for x in xx], '-', lw=2, label='bspline-base-1')
ax.plot(xx, [calc_bspline(x, 2, N, M, x_max, x_min) * c[2] for x in xx], '-', lw=2, label='bspline-base-2')
ax.plot(x_min - h / 2.0, c[0], '.')
ax.plot(x_min + 1.0 * h / 2.0, c[1], '.')
ax.plot(x_min + 3.0 * h / 2.0, c[2], '.')
ax.legend(loc='best')
plt.show()
"""
