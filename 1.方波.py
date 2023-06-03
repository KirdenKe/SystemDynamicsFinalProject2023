import numpy as np
import matplotlib.pylab as plt
from scipy.integrate import odeint

Posmag = 2.5
Negmag = -2.5
amplitude = Posmag - Negmag
ncycle = 1
frequence = 0.100234
period = 1 / frequence
totaltime = ncycle * period

dt = 0.001
npoint = totaltime / dt

def ODE(x, t, m1, m2, b1, b2, k1, k2, F):
    # 方波
    a = t / period
    aint = int(a)
    residual = a - aint
    if (residual <= 0.5):
        F = Posmag
        dx = x[1]; dx2 = x[2]; dx3 = x[3]
        dx4 = (k2 * F - (m2 * (b1 + b2) + m1 * b2) * x[3] - (m2 * (k1 + k2) + b1 * b2 + m1 * k2) * x[2] - (b1 * k2 + b2 * k1) * x[1] - k1 * k2 * x[0]) / (m1 * m2)
        return [dx, dx2, dx3, dx4]
    else:
        F = Negmag
        dx = x[1]; dx2 = x[2]; dx3 = x[3]
        dx4 = (k2 * F - (m2 * (b1 + b2) + m1 * b2) * x[3] - (m2 * (k1 + k2) + b1 * b2 + m1 * k2) * x[2] - (b1 * k2 + b2 * k1) * x[1] - k1 * k2 * x[0]) / (m1 * m2)
        return [dx, dx2, dx3, dx4]

t = np.arange(0, totaltime + dt, dt)
square = np.zeros(t.size)
time = 0
for i in range(0, t.size):
    a = time / period
    a_int = int(a)
    residual = a - a_int
    if (residual <= 0.5):
        square[i] = Posmag
    else:
        square[i] = Negmag
    time += dt

m1 = 0.1
m2 = 0.001
b1 = 0.1
b2 = 0.001
k1 = 1
k2 = 1

x_initial = [0, 0, 0, 0]
x = odeint(ODE, x_initial, t, args = (m1, m2, b1, b2, k1, k2, square))
plt.plot(t, square)
plt.plot(t, x[:, 0])
plt.grid()
plt.show()