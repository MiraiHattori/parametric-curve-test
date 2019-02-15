import numpy as np
import matplotlib.pyplot as plt
import copy
import random

# 一般のB-Spline
# this list has to be long enough
#colorlist = ["red", "green", "blue", "yellow", "cornflowerblue", "darkorange", "magenta", "black"]
colorlist_base_before = ["paleturquoise"] * 10
colorlist_base_after = ["darksalmon"] * 10

colorlist_sum_before = ["darkturquoise"] * 10
colorlist_sum_after = ["orangered"] * 10

colorlist_point_before = ["darkturquoise"] * 10
colorlist_point_after = ["orangered"] * 10

colorlist_convex_hull_before = ["darkturquoise"] * 10
colorlist_convex_hull_after = ["orangered"] * 10


def calc_bspline(t, i, n, t_knot):
    if len(t_knot) <= i+n+1:
        print("Error: len(t_knot) (= {} ) is invalid".format(len(t_knot)))
        return 0
    return calc_bspline_base(t, i, n, n, t_knot)

def calc_bspline_base(t, i, n, _n, t_knot):
    ti = t_knot[i]
    ti_1 = t_knot[i+1]
    ti_n = t_knot[i+_n]
    ti_n_1 = t_knot[i+_n+1]
    if _n == 0:
        if ti <= t and t <= ti_1:
            return 1
        else:
            return 0
    else:
        return ((t - ti) / (ti_n - ti) * calc_bspline_base(t, i, n, _n - 1, t_knot) + \
               (ti_n_1 - t) / (ti_n_1 - ti_1) * calc_bspline_base(t, i + 1, n, _n - 1, t_knot))

n = 3
t_knot = [0, 2, 3, 6, 9, 10]
m = len(t_knot)
c = [1.0] * (m - n - 1)
fig, ax = plt.subplots(num=0)
tt = np.linspace(t_knot[0], t_knot[-1], 999)

for i in range(len(c)):
    ax.plot(tt, [calc_bspline(t, i, n, t_knot) * c[i] for t in tt], '-', lw=2, label='bspline-base-'+str(i))
sum_bsp = []
for t in tt:
    s = 0.0
    for i in range(m - n - 1):
        s += c[i] * calc_bspline(t, i, n, t_knot)
    sum_bsp.append(s)
ax.plot(tt, sum_bsp, '-', lw=2, label='bspline-sum')
t_used = []
t_unused = []
for i in range(0, m):
    if n <= i and i <= m - n - 1:
        t_used.append(t_knot[i])
    else:
        t_unused.append(t_knot[i])
ax.vlines(t_used, 0, max(c), "red", linestyles='dashed')
ax.vlines(t_unused, 0, max(c), "blue", linestyles='dashed')
ax.legend(loc='best')
#plt.show()



# 一様B-Spline
def calc_cardinal_bspline(t, i, n, m, t_max, t_min):
    return calc_cardinal_bspline_base(t, i, n, n, m, t_max, t_min)

def calc_cardinal_bspline_base(t, i, n, _n, m, t_max, t_min):
    h = (t_max - t_min) * 1.0 / (m - n)
    ti = h * i + (m * t_min - n * t_max) * 1.0 / (m - n)
    ti_1 = ti + h
    ti_n_1 = ti_1 + _n * h
    if _n == 0:
        if ti <= t and t <= ti_1:
            return 1
        else:
            return 0
    else:
        return ((t - ti) * calc_cardinal_bspline_base(t, i, n, _n - 1, m, t_max, t_min) + \
               (ti_n_1 - t) * calc_cardinal_bspline_base(t, i + 1, n, _n - 1, m, t_max, t_min)) / (_n * h)


n = 3
m = 6
t_max = 1
t_min = 0

h = (t_max - t_min) * 1.0 / (m - n)
c = [1.0] * m

fig, ax = plt.subplots(num=1)
tt = np.linspace((m * t_min - n * t_max) * 1.0 / (m - n), h * (m + n) + (m * t_min - n * t_max) * 1.0 / (m - n), 999)

for i in range(m):
    ax.plot(tt, [calc_cardinal_bspline(t, i, n, m, t_max, t_min) * c[i] for t in tt], '-', lw=2, label='bspline-base-'+str(i))

sum_bsp = []
for t in tt:
    s = 0.0
    for i in range(m):
        s += c[i] * calc_cardinal_bspline(t, i, n, m, t_max, t_min)
    sum_bsp.append(s)
ax.plot(tt, sum_bsp, '-', lw=2, label='bspline-sum')

t_used = []
t_unused = []

for i in range(0, m + n + 1):
    if n <= i and i <= m:
        t_used.append(h * i + (m * t_min - n * t_max) / (m - n))
    else:
        t_unused.append(h * i + (m * t_min - n * t_max) / (m - n))
ax.vlines(t_used, 0, max(c), "red", linestyles='dashed')
ax.vlines(t_unused, 0, max(c), "blue", linestyles='dashed')
for i in range(m):
    ax.plot(h * (i + i + n + 1) / 2.0 + (m * t_min - n * t_max) / (m - n), c[i], '.')
#ax.legend(loc='best')
#plt.show()

# 局所性
n = 3
m = 6
t_max = 1
t_min = 0

h = (t_max - t_min) * 1.0 / (m - n)
c = [1.0] * m
# i=n+1番目を変更
c_modified = copy.deepcopy(c)
c_modified[n+1] = 1.4

fig, ax = plt.subplots(num=2)
tt = np.linspace((m * t_min - n * t_max) * 1.0 / (m - n), h * (m + n) + (m * t_min - n * t_max) * 1.0 / (m - n), 999)

for i in range(m):
    ax.plot(tt, [calc_cardinal_bspline(t, i, n, m, t_max, t_min) * c[i] for t in tt], '-', lw=2, label='bspline-base-'+str(i), color=colorlist_base_before[i])

for i in range(m):
    ax.plot(tt, [calc_cardinal_bspline(t, i, n, m, t_max, t_min) * c_modified[i] for t in tt], '-', lw=2, label='bspline-base-modified-'+str(i), color=colorlist_base_after[i])

sum_bsp = []
for t in tt:
    s = 0.0
    for i in range(m):
        s += c[i] * calc_cardinal_bspline(t, i, n, m, t_max, t_min)
    sum_bsp.append(s)
ax.plot(tt, sum_bsp, '-', lw=2, label='bspline-sum', color=colorlist_sum_before[m+1])

sum_bsp = []
for t in tt:
    s = 0.0
    for i in range(m):
        s += c_modified[i] * calc_cardinal_bspline(t, i, n, m, t_max, t_min)
    sum_bsp.append(s)
ax.plot(tt, sum_bsp, '-', lw=2, label='bspline-sum-modified', color=colorlist_sum_after[m+1])
t_used = []
t_unused = []

for i in range(0, m + n + 1):
    if n <= i and i <= m:
        t_used.append(h * i + (m * t_min - n * t_max) / (m - n))
    else:
        t_unused.append(h * i + (m * t_min - n * t_max) / (m - n))
ax.vlines(t_used, 0, max(c + c_modified) * 1.5, "red", linestyles='dashed')
ax.vlines(t_unused, 0, max(c + c_modified) * 1.5, "blue", linestyles='dashed')

for i in range(m):
    ax.plot(h * (i + i + n + 1) / 2.0 + (m * t_min - n * t_max) / (m - n), c[i], '.', markersize=10, color=colorlist_point_before[i])

for i in range(m):
    ax.plot(h * (i + i + n + 1) / 2.0 + (m * t_min - n * t_max) / (m - n), c_modified[i], '.', markersize=10, color=colorlist_point_after[i])
#ax.legend(loc='best')
#plt.show()


# 凸包性
n = 3
m = 6
t_max = 1
t_min = 0

h = (t_max - t_min) * 1.0 / (m - n)

fig, ax = plt.subplots(num=3)
tt = np.linspace((m * t_min - n * t_max) * 1.0 / (m - n), h * (m + n) + (m * t_min - n * t_max) * 1.0 / (m - n), 999)

for i in range(m):
    c_modified[i] = random.uniform(0.0, 1.5)

for i in range(m):
    ax.plot(tt, [calc_cardinal_bspline(t, i, n, m, t_max, t_min) * c_modified[i] for t in tt], '-', lw=2, label='bspline-base-modified-'+str(i), color=colorlist_base_after[i])

sum_bsp = []
for t in tt:
    s = 0.0
    for i in range(m):
        s += c_modified[i] * calc_cardinal_bspline(t, i, n, m, t_max, t_min)
    sum_bsp.append(s)
ax.plot(tt, sum_bsp, '-', lw=2, label='bspline-sum-modified', color=colorlist_sum_after[m+1])

c_modified_minus = []
c_modified_plus = []

for i in range(m):
    t = h * (n + i - 1) + (m * t_min - n * t_max) / (m - n)
    s = 0.0
    for j in range(m):
        s += c_modified[j] * calc_cardinal_bspline(t+0.001, j, n, m, t_max, t_min)
    if c_modified[i] >= s:
        c_modified_plus.append([t, c_modified[i]])
    else:
        c_modified_minus.append([t, c_modified[i]])
ax.plot([cp[0] for cp in c_modified_plus], [cp[1] for cp in c_modified_plus], '-', lw=2, label='bspline-convex-hull-upper', color=colorlist_convex_hull_after[m+1])
ax.plot([cm[0] for cm in c_modified_minus], [cm[1] for cm in c_modified_minus], '-', lw=2, label='bspline-convex-hull-lower', color=colorlist_convex_hull_after[m+1])
if len(c_modified_minus) > 0 and len(c_modified_plus) > 0:
    ax.plot([c_modified_plus[0][0], c_modified_minus[0][0]], [c_modified_plus[0][1], c_modified_minus[0][1]], '-', lw=2, label='bspline-convex-hull-lower-upper', color=colorlist_convex_hull_after[m+1])
    ax.plot([c_modified_plus[-1][0], c_modified_minus[-1][0]], [c_modified_plus[-1][1], c_modified_minus[-1][1]], '-', lw=2, label='bspline-convex-hull-upper-lower', color=colorlist_convex_hull_after[m+1])

t_used = []
t_unused = []

for i in range(0, m + n + 1):
    if n <= i and i <= m:
        t_used.append(h * i + (m * t_min - n * t_max) / (m - n))
    else:
        t_unused.append(h * i + (m * t_min - n * t_max) / (m - n))
ax.vlines(t_used, 0, max(c + c_modified) * 1.5, "red", linestyles='dashed')
ax.vlines(t_unused, 0, max(c + c_modified) * 1.5, "blue", linestyles='dashed')

#for i in range(m):
#    ax.plot(h * (i + i + n + 1) / 2.0 + (m * t_min - n * t_max) / (m - n), c[i], '.', markersize=10, color=colorlist_point_before[i])

for i in range(m):
    ax.plot(h * (i + i + n + 1) / 2.0 + (m * t_min - n * t_max) / (m - n), c_modified[i], '.', markersize=10, color=colorlist_point_after[i])
#ax.legend(loc='best')
plt.show()
