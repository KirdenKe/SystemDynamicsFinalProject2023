import numpy as np
import matplotlib.pylab as plt
from scipy import signal
import math

m1 = 1
m2 = 10
b1 = 1
b2 = 10
k1 = 1
k2 = 1

num = [k2]
den = [m1 * m2, m2 * (b1 + b2) + m1 * b2, m2 * (k1 + k2) + b1 * b2 + m1 * k2, b1 * k2 + b2 * k1, k1 * k2]
system = signal.lti(num, den)
t, y = signal.step(system)
plt.plot(t, y)
plt.grid()
plt.show()

f = np.logspace(-2, 2, num=1000)  # 頻率
w = 2 * np.pi * f  # 角頻率
w, mag, phase = signal.bode(system, w)
for i in range(len(phase) - 1):
    if (phase[i] + 90) * (phase[i + 1] + 90) < 0:
        print("Resonance: {0:.5f}".format((phase[i + 1] + 90) / (phase[i + 1] - phase[i]) * (f[i + 1] - f[i]) + f[i])) # 共振頻率
    elif (phase[i] + 180) * (phase[i + 1] + 180) < 0:
        print("{0:.5f}".format((phase[i + 1] + 180) / (phase[i + 1] - phase[i]) * (f[i + 1] - f[i]) + f[i]))

plt.subplot(211)
plt.semilogx(f, mag)
plt.xlabel('frequency(1/sec)')
plt.ylabel('magnitude(db)')
plt.grid()

plt.subplot(212)
plt.semilogx(f, phase)
plt.xlabel('frequency(1/sec)')
plt.ylabel('phase(degree)')
plt.grid()
plt.show()