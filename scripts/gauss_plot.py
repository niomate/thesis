import numpy as np
import matplotlib.pyplot as plt

epsilon = 0.001


@np.vectorize
def gauss(x, sigma):
    return 1/(sigma * np.sqrt(2*np.pi)) * np.exp(-x**2 / (2*sigma**2))


def end(sigma, prec=5):
    return int(np.floor(sigma * prec)) + 1


def xarr(sigma, prec=5):
    return np.array(list(np.arange(-end(sigma, prec), end(sigma, prec)+1, 1)))


def plot_gauss(sigma, prec=5):
    x = xarr(sigma, prec)
    y = gauss(x, sigma)
    print(np.sum(y))
    line1 = np.argmin(y < epsilon)
    line2 = np.argmax((x > 0) & (y < epsilon))
    plt.plot(x, y)
    plt.axvline(x=sigma, c='b')
    plt.axvline(x=-sigma, c='b')
    plt.axvline(x=x[line1], c='r')
    plt.axvline(x=x[line2], c='r')
    plt.grid(b=True)
    plt.show()


def calc_ratio(sigma):
    x = xarr(sigma, 5)
    y = gauss(x, sigma)
    line = np.argmax((x > 0) & (y < epsilon))
    print("Sigma:", sigma, "first y < epsilon:", x[line], "Ratio:", x[line] / sigma)
    return x[line] / sigma

# sigmas = np.arange(1, 10, 0.5)
# ratios = []
# for sigma in sigmas:
    # r = calc_ratio(sigma)
    # if np.isclose(r, 0):
        # break
    # ratios.append(r)

# plt.plot(ratios, sigmas[:len(ratios)])
# plt.show()
