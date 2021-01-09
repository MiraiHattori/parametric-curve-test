import numpy as np
import matplotlib.pyplot as plt
import math

fig, ax = plt.subplots(num=3)

def mod2sgn(n):
    return 1 if n % 2 == 0 else -1


def factorial(n):
    return 1 if n < 2 else n * factorial(n - 1)

# TODO trapozoidal integral
def fresnel_sin(t_start, t_end, n):
    t = np.linspace(t_start, t_end, n)
    if n == 1:
        print("in fresnel_sin: n is invalid {}".format(n))
        return 0.0
    return np.sum(np.sin(t * t)) * (t_end - t_start) / (n - 1)


def fresnel_cos(u, n):
    ans = 0.0
    for i in range(n):
        ans += mod2sgn(i) * math.pow(u, 4 * i + 1) / ((4 * i + 1) * factorial(2 * i))
    return ans


print(fresnel_sin(0.0, 5, 1000))
#
#u = np.linspace(0, 6, 1000)
#t = np.vectorize(lambda u: fresnel_cos(0.0, u, 3, 1000))(u)
#x = np.vectorize(lambda u: fresnel_sin(0.0, u, 3, 1000))(u)
#
#ax.plot(t, x, '.', markersize=10)
##ax.legend(loc='best')
#plt.show()
