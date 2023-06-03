import numpy as np
import matplotlib.pylab as plt
from scipy import signal
import math

def DampingRatio(x, tar): # 參考至https://electronics.stackexchange.com/questions/565679/finding-resonant-frequency-or-damping-ratio-from-bode-plot
    return np.abs((np.sqrt(1 - 2 * np.power(x, 2) + np.sqrt(2 - 4 * np.power(x, 2) + 4 * np.power(x, 4)))) - tar)

m1 = 1.6
m2 = 1
b1 = 1
b2 = 1
k1 = 1
k2 = 1

num = [b2, k2]
den = [m1 * m2, m2 * (b1 + b2) + m1 * b2, m2 * (k1 + k2) + b1 * b2 + m1 * k2, b1 * k2 + b2 * k1, k1 * k2]
system = signal.lti(num, den)
t, y = signal.step(system)
plt.plot(t, y)
plt.grid()
plt.show()

f = np.logspace(-2, 2, num = 1000)  # 頻率
w = 2 * np.pi * f  # 角頻率
w, mag, phase = signal.bode(system, w)
fn = 0; MagN = 0
for i in range(len(phase) - 1):
    if (phase[i] + 90) * (phase[i + 1] + 90) < 0:
        print("Natural: {0:.5f}".format((phase[i + 1] + 90) / (phase[i + 1] - phase[i]) * (f[i + 1] - f[i]) + f[i])) # 自然頻率
    elif (phase[i] + 180) * (phase[i + 1] + 180) < 0:
        MagN = (phase[i + 1] + 180) / (phase[i + 1] - phase[i]) * (mag[i + 1] - mag[i]) + mag[i]
        fn = (phase[i + 1] + 180) / (phase[i + 1] - phase[i]) * (f[i + 1] - f[i]) + f[i]
        print("Resonance: {0:.5f}".format(fn))
        break

for i in range(len(mag) - 1):
    if (mag[i] - MagN + 3) * (mag[i + 1] - MagN + 3) < 0:
        f3db = (mag[i + 1] - MagN + 3) / (mag[i + 1] - mag[i]) * (f[i + 1] - f[i]) + f[i]
        if f3db > fn:
            break

gr = 0.5 * (3 - np.sqrt(5)) # golden ratio
eps = 1e-5 # allowable error
xL = 0; xR = 2.5
x1 = (1 - gr) * xL + gr * xR; x2 = gr * xL + (1 - gr) * xR;
f1 = DampingRatio(x1, f3db / fn); f2 = DampingRatio(x2, f3db / fn)
while xR - xL > eps:
    if f1 < f2:
        xR = x2; x2 = x1; f2 = f1
        x1 = (1 - gr) * xL + gr * xR
        f1 = DampingRatio(x1, f3db / fn)
    else:
        xL = x1; x1 = x2; f1 = f2;
        x2 = gr * xL + (1 - gr) * xR
        f2 = DampingRatio(x2, f3db / fn)
print("Damping Ratio: {0:.3f}".format(0.5 * (xL + xR)))

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