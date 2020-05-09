import matplotlib.pyplot as plt
import numpy as np


def weickert(s, l=1.0):
    if s**2 == 0:
        return 1
    bot = (s / l) ** 8
    v = 1 - np.exp(-3.31488 / bot)
    return v


weickert = np.vectorize(weickert, otypes=[float])


def charbonnier(s, l=1.0):
    return 1 / np.sqrt(1 + (s/l)**2)


charbonnier = np.vectorize(charbonnier, otypes=[float])

x = np.arange(0, 5, 0.01)
plt.plot(x, weickert(x), 'g')
plt.plot(x, charbonnier(x), 'r')
plt.grid(True)
plt.axvline(x=1, linestyle='--')
plt.legend(['Weickert', 'Charbonnier', 'Lambda'])
plt.show()
