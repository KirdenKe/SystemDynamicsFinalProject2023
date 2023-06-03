import numpy as np
import matplotlib.pylab as plt
from scipy.integrate import odeint
from vpython import *

Posmag = 2.5
Negmag = -2.5
amplitude = Posmag - Negmag
ncycle = 5
frequence = 0.17216
period = 1 / frequence
totaltime = ncycle * period

dt = 0.001
npoint = totaltime / dt

def ODE(x, t, m1, m2, b1, b2, k1, k2, F):
    # 簡諧波
    a = t / period
    aint = int(a)
    residual = a - aint
    F = Posmag * np.sin(residual * 2 * np.pi)
    dF = Posmag * (2 * np.pi * frequence) * np.cos(residual * 2 * np.pi)
    dx = x[1];
    dx2 = x[2];
    dx3 = x[3]
    dx4 = (b2 * dF + k2 * F - (m2 * (b1 + b2) + m1 * b2) * x[3] - (m2 * (k1 + k2) + b1 * b2 + m1 * k2) * x[2] - (b1 * k2 + b2 * k1) * x[1] - k1 * k2 * x[0]) / (m1 * m2)
    return [dx, dx2, dx3, dx4]

t = np.arange(0, totaltime + dt, dt)
harmonic = np.zeros(t.size)
time = 0
for i in range(0, t.size):
    a = time / period
    a_int = int(a)
    residual = a - a_int
    harmonic[i] = Posmag * np.sin(residual * 2 * np.pi)
    time += dt

m1 = 1.6
m2 = 1
b1 = 1
b2 = 1
k1 = 1
k2 = 1

x_initial = [0, 0, 0, 0]
x = odeint(ODE, x_initial, t, args = (m1, m2, b1, b2, k1, k2, harmonic))
plt.plot(t, harmonic)
plt.plot(t, x[:, 0])
plt.grid()
plt.show()