import numpy as np
import matplotlib.pylab as plt
from scipy.integrate import odeint
from vpython import *

Posmag = 2.5
Negmag = -2.5
amplitude = Posmag - Negmag
ncycle = 5
frequence = 0.02874
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

def ODE2(x, t, m1, m2, b1, b2, k1, k2, F):
    # 方波
    a = t / period
    aint = int(a)
    residual = a - aint
    if (residual <= 0.5):
        F = Posmag
        dx = x[1]; dx2 = x[2]; dx3 = x[3]
        dx4 = ((k1 + k2) * F - (m2 * (b1 + b2) + m1 * b2) * x[3] - (m2 * (k1 + k2) + b1 * b2 + m1 * k2) * x[2] - (b1 * k2 + b2 * k1) * x[1] - k1 * k2 * x[0]) / (m1 * m2)
        return [dx, dx2, dx3, dx4]
    else:
        F = Negmag
        dx = x[1]; dx2 = x[2]; dx3 = x[3]
        dx4 = ((k1 + k2) * F - (m2 * (b1 + b2) + m1 * b2) * x[3] - (m2 * (k1 + k2) + b1 * b2 + m1 * k2) * x[2] - (b1 * k2 + b2 * k1) * x[1] - k1 * k2 * x[0]) / (m1 * m2)
        return [dx, dx2, dx3, dx4]

L = 0.5
scene = canvas(width=800, height=500, center=vector(0, -L * 0.8, 0), range=1.2 * L)  # 設定畫面
ceiling = box(length=1.2, height=0.005, width=0.4, opacity=1, color=color.blue, shininess=1)  # 畫天花板

floor = box(length=1.2, height=0.005, width=0.4, opacity=0.2)  # 畫地板
floor.pos = vector(0, -L*1.5, 0)

wall = box(length=0.005, height=0.4, width=0.4, color=color.green, opacity=0.5)  # 畫牆
wall.pos = vector(-0.5,-L*1.1, 0)

object1 = box(length=0.15, height=0.15, width=0.15, color=color.white, opacity=1)  # 畫m1
object1.pos = vector(-0.025, -1.35*L, 0)

spring1 = helix(radius=0.02, thickness=0.01)  # 畫k1
spring1.pos = vector(-0.5, -1.28*L, 0)  # 彈簧頭端的位置
spring1.color = vector(0.7, 0.5, 0.2)
spring1.axis = 0.4 * vector(1,0,0)

damperB1 = cylinder(radius=0.02, thickness=0.01)
damperB1.pos = vector(-0.5, -1.45*L, 0)
damperB1.color=vector(0.3, 0.1, 0.3)
damperB1.axis = 0.2 * vector(1, 0, 0)

damperM1 = cylinder(radius=0.01, thickness=0.01)
damperM1.pos = vector(-0.5, -1.45*L, 0)
damperM1.color = vector(0.1, 0.7, 0.1)
damperM1.axis = 0.4 * vector(1, 0, 0)

object2 = box(length=0.15, height=0.15, width=0.15, color=color.white, opacity=1)  # 畫m2
object2.pos = vector(0.525, -1.35*L, 0)

spring2 = helix(radius=0.02, thickness=0.01)  # 畫k2
spring2.pos = vector(0.05, -1.28*L, 0)  # 彈簧頭端的位置
spring2.color = vector(0.7, 0.5, 0.2)
spring2.axis = 0.4 * vector(1,0,0)

damperB2 = cylinder(radius=0.02, thickness=0.01)
damperB2.pos = vector(0.05, -1.45*L, 0)
damperB2.color=vector(0.3, 0.1, 0.3)
damperB2.axis = 0.2 * vector(1, 0, 0)

damperM2 = cylinder(radius=0.01, thickness=0.01)
damperM2.pos = vector(0.05, -1.45*L, 0)
damperM2.color = vector(0.1, 0.7, 0.1)
damperM2.axis = 0.4 * vector(1, 0, 0)

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

m1 = 1
m2 = 10
b1 = 1
b2 = 10
k1 = 1
k2 = 1

x_initial = [0, 0, 0, 0]
x1 = odeint(ODE, x_initial, t, args = (m1, m2, b1, b2, k1, k2, square))

x_initial = [0, 0, 0, 0]
x2 = odeint(ODE2, x_initial, t, args = (m1, m2, b1, b2, k1, k2, square))

j = 0

while True:
    rate(1 / dt)
    object1.pos = vector(-0.025 + x1[j, 0] / 50, -1.35 * L, 0)
    spring1.axis = vector(object1.pos.x + 0.5 - 0.075, 0, 0)
    damperM1.axis = vector(object1.pos.x + 0.5 - 0.075, 0, 0)
    spring2.pos = vector(object1.pos.x + 0.075, -1.28*L, 0)
    damperB2.pos = vector(object1.pos.x + 0.075, -1.45 * L, 0)
    damperM2.pos = vector(object1.pos.x + 0.075, -1.45 * L, 0)
    object2.pos = vector(0.525 + x2[j, 0] / 50, -1.35 * L, 0)
    spring2.axis = vector(object2.pos.x - object1.pos.x - 0.15, 0, 0)
    damperM2.axis = vector(object2.pos.x - object1.pos.x - 0.15, 0, 0)
    j += 1
    test = x1[:, 0].size
    if j > x1[:, 0].size - 1:
        j = 0