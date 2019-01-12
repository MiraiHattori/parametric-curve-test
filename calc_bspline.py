import numpy as np
import matplotlib.pyplot as plt

def calc_bspline(t, i, n, t_knot):
    if len(t_knot) <= i+n+1:
        print("Error: len(t_knot) (= {} ) is invalid".format(len(t_knot)))
        return 0
    return calc_bspline_base(t, i, n, n, t_knot)

def calc_bspline_base(t, i, n, _n, t_knot):
    ti = t_knot[i]
    ti_1 = t_knot[i+1]
    ti_n_1 = t_knot[i+n+1]
    if t < ti or ti_1 <= t:
        return 0
    elif _n == 0:
        return 1
    else:
        return ((t - ti) * calc_cardinal_bspline_base(t, i, n, _n - 1, t_knot) + \
               (ti_n_1 - t) * calc_cardinal_bspline_base(t, i + 1, n, _n - 1, t_knot)) / (_n * h)

n = 3
t_knot = [1, 2, 3, 4, 5, 6]
c = [1.0] * len(t_knot)
#
#
#fig, ax = plt.subplots()
#tt = np.linspace(t_min - h * 2, t_max + h * 2, 999)
#
#t = []
#
#for i in range(m):
#    ax.plot(tt, [calc_cardinal_bspline(t, i, n, m, t_max, t_min) * c[i] for t in tt], '-', lw=2, label='bspline-base-'+str(i))
#
#sum_bsp = []
#for t in tt:
#    s = 0.0
#    for i in range(m):
#        s += calc_cardinal_bspline(t, i, n, m, t_max, t_min)
#    sum_bsp.append(s)
#ax.plot(tt, sum_bsp, '-', lw=2, label='bspline-base-'+str(i))
#ax.legend(loc='best')
#plt.show()



#def calc_cardinal_bspline(t, i, n, m, t_max, t_min):
#    return calc_cardinal_bspline_base(t, i, n, n, m, t_max, t_min)
#
#def calc_cardinal_bspline_base(t, i, n, _n, m, t_max, t_min):
#    h = (t_max - t_min) * 1.0 / (m - n)
#    ti = h * i + (m * t_min - n * t_max) * 1.0 / (m - n)
#    ti_1 = ti + h
#    ti_n_1 = ti_1 + _n * h
#    if _n == 0:
#        if t < ti or ti_1 < t:
#            return 0
#        elif _n == 0:
#            return 1
#    else:
#        return ((t - ti) * calc_cardinal_bspline_base(t, i, n, _n - 1, m, t_max, t_min) + \
#               (ti_n_1 - t) * calc_cardinal_bspline_base(t, i + 1, n, _n - 1, m, t_max, t_min)) / (_n * h)
#
#
#n = 3
#m = 14
#t_max = 1
#t_min = 0
#
#h = (t_max - t_min) * 1.0 / (m - n)
#c = [1.0] * m
#
#fig, ax = plt.subplots()
#tt = np.linspace(, t_max + h * 2, 999)
#
#t = []
#
#for i in range(m):
#    ax.plot(tt, [calc_cardinal_bspline(t, i, n, m, t_max, t_min) * c[i] for t in tt], '-', lw=2, label='bspline-base-'+str(i))
#
#sum_bsp = []
#for t in tt:
#    s = 0.0
#    for i in range(m):
#        s += calc_cardinal_bspline(t, i, n, m, t_max, t_min)
#    sum_bsp.append(s)
#ax.plot(tt, sum_bsp, '-', lw=2, label='bspline-base-'+str(i))
#ax.legend(loc='best')
#plt.show()
