import numpy as np
import matplotlib.pyplot as plt

def B(x, k, i, t):
   if k == 0:
      return 1.0 if t[i] <= x < t[i+1] else 0.0
   if t[i+k] == t[i]:
      c1 = 0.0
   else:
      c1 = (x - t[i])/(t[i+k] - t[i]) * B(x, k-1, i, t)
   if t[i+k+1] == t[i+1]:
      c2 = 0.0
   else:
      c2 = (t[i+k+1] - x)/(t[i+k+1] - t[i+1]) * B(x, k-1, i+1, t)
   return c1 + c2

def bspline(x, t, c, k):
   n = len(t) - k - 1
   assert (n >= k+1) and (len(c) >= n)
   return sum(c[i] * B(x, k, i, t) for i in range(n))

k = 2
t = [-1, 0, 1, 2, 3, 4, 5]
c = [2.1, 2.2, 3.2, 1.2, 2.2]
c_modified = [2.1, 2.2, 4.2, 1.2, 2.2]
fig, ax = plt.subplots()
xx = np.linspace(-1, 8, 100)
# ax.plot(xx, [bspline(x, t, c, k) for x in xx], 'r-', lw=2, label='bspline')
ax.plot(xx, [B(x, k, 0, t) * c[0] for x in xx], '-', lw=2, label='bspline-base-0')
ax.plot(xx, [B(x, k, 1, t) * c[1] for x in xx], '-', lw=2, label='bspline-base-1')
ax.plot(xx, [B(x, k, 2, t) * c[2] for x in xx], '-', lw=2, label='bspline-base-2')
ax.plot(xx, [B(x, k, 3, t) * c[3] for x in xx], '-', lw=2, label='bspline-base-3')
ax.plot(xx, [B(x, k, 0, t) * c[0] + B(x, k, 1, t) * c[1] + B(x, k, 2, t) * c[2] + B(x, k, 3, t) * c[3]  for x in xx], '-', lw=2, label='bspline-base-sum', linestyle="dashed")
ax.plot(xx, [B(x, k, 0, t) * c[0] + B(x, k, 1, t) * c[1] + B(x, k, 2, t) * c_modified[2] + B(x, k, 3, t) * c[3]  for x in xx], '-', lw=2, label='bspline-base-sum-modified')
ax.plot(0.5, c[0], '.')
ax.plot(1.5, c[1], '.')
ax.plot(2.5, c[2], '.')
ax.plot(2.5, c_modified[2], '.')
ax.plot(3.5, c[3], '.')
# ax.plot(4.5, c[4], '.')
ax.legend(loc='best')
plt.show()
