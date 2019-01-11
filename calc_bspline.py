import numpy as np
import matplotlib.pyplot as plt

def calc_bspline(x, j, N, M, x_max, x_min):
    return calc_bspline_base(x, j, N, N, M, x_max, x_min)

def calc_bspline_base(x, j, N, _n, M, x_max, x_min):
    h = (x_max - x_min) * 1.0 / (M - N)
    xj = h * j + (M * x_min - N * x_max) * 1.0 / (M - N)
    xj_1 = xj + h
    xj_n_1 = xj_1 + _n * h
    if x < xj or xj_n_1 <= x:
        return 0
    elif _n == 0:
        return 1
    else:
        return ((x - xj) * calc_bspline_base(x, j, N, _n - 1, M, x_max, x_min) + \
               (xj_n_1 - x) * calc_bspline_base(x, j + 1, N, _n - 1, M, x_max, x_min)) / (_n * h)


N = 5
M = 14
x_max = 1
x_min = 0

h = (x_max - x_min) * 1.0 / (M - N)
c = []
for i in range(M):
    c.append(i * 1.0 + 1)

t_used = []
t_unused = []

for i in range(0, M + N + 1):
    if N <= i and i <= M:
        t_used.append(h * i + (M * x_min - N * x_max) / (M - N))
    else:
        t_unused.append(h * i + (M * x_min - N * x_max) / (M - N))

fig, ax = plt.subplots()
xx = np.linspace(h * 0 + (M * x_min - N * x_max) / (M - N), h * (M + N) + (M * x_min - N * x_max) / (M - N), 999)

for i in range(M):
    ax.plot(xx, [calc_bspline(x, i, N, M, x_max, x_min) * c[i] for x in xx], '-', lw=2, label='bspline-base-'+str(i))

sum_bsp = []
for x in xx:
    s = 0.0
    for i in range(M):
        s += c[i] * calc_bspline(x, i, N, M, x_max, x_min)
    sum_bsp.append(s)
ax.plot(xx, sum_bsp, '-', lw=2, label='bspline-sum')
ax.vlines(t_used, 0, max(c), "red", linestyles='dashed')
ax.vlines(t_unused, 0, max(c), "blue", linestyles='dashed')
for i in range(M):
    ax.plot(h * (i + i + N + 1) / 2.0 + (M * x_min - N * x_max) / (M - N), c[i], '.')
#ax.plot(x_min - h / 2.0, c[0], '.')
#ax.plot(x_min + 1.0 * h / 2.0, c[1], '.')
#ax.plot(x_min + 3.0 * h / 2.0, c[2], '.')
ax.legend(loc='best')
plt.show()
