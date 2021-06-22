import numpy as np

from scipy import integrate
from matplotlib import pyplot as plt

mu1 = np.array([-2.0, -2.0])
mu2 = np.array([1.0, 1.0])

std1 = np.array([1.5, 0.5])
std2 = np.array([0.5, 1.5])


def gaussian(X, Y, mu, std):
    assert len(mu) == len(std) == 2
    return np.exp(-((X - mu[0]) ** 2) / std[0] ** 2 - (Y - mu[1]) ** 2 / std[1] ** 2)


delta = 0.1
x = np.arange(-5.0, 5.0, delta)
y = np.arange(-5.0, 5.0, delta)

X, Y = np.meshgrid(x, y)

d1 = gaussian(X, Y, mu1, std1)
d2 = gaussian(X, Y, mu2, std2)

fig = plt.figure()
ax = plt.subplot(1, 2, 1)
ax.set(aspect="equal")
ax.contour(X, Y, d1)
ax = plt.subplot(1, 2, 2)
ax.set(aspect="equal")
ax.contour(X, Y, d2)
plt.savefig("images/contours.png")


def fx(x, y, mu, std):
    return x * gaussian(x, y, mu, std)


def fy(x, y, mu, std):
    return y * gaussian(x, y, mu, std)


def fxy_fyx(x, y, mu, std):
    return x * y * gaussian(x, y, mu, std)


def fxx(x, y, mu, std):
    return x ** 2 * gaussian(x, y, mu, std)


def fyy(x, y, mu, std):
    return y ** 2 * gaussian(x, y, mu, std)


def density_integral(x, y, mu, std):
    I = integrate.dblquad(
        lambda x, y: gaussian(x, y, mu, std),
        -np.inf,
        np.inf,
        lambda x: -np.inf,
        lambda x: np.inf,
    )

    return I[0]


def shape_centroid(mu, std):
    Rx = integrate.dblquad(
        lambda x, y: fx(x, y, mu, std),
        -np.inf,
        np.inf,
        lambda x: -np.inf,
        lambda x: np.inf,
    )

    Ry = integrate.dblquad(
        lambda x, y: fy(x, y, mu, std),
        -np.inf,
        np.inf,
        lambda x: -np.inf,
        lambda x: np.inf,
    )

    # Denominator
    R0 = density_integral(x, y, mu, std)

    return np.array([Rx[0], Ry[0]]) / R0


def quadrupole(mu, std):
    Mxx = integrate.dblquad(
        lambda x, y: fxx(x, y, mu, std),
        -np.inf,
        np.inf,
        lambda x: -np.inf,
        lambda x: np.inf,
    )

    Mxy = Myx = integrate.dblquad(
        lambda x, y: fxy_fyx(x, y, mu, std),
        -np.inf,
        np.inf,
        lambda x: -np.inf,
        lambda x: np.inf,
    )

    Myy = integrate.dblquad(
        lambda x, y: fyy(x, y, mu, std),
        -np.inf,
        np.inf,
        lambda x: -np.inf,
        lambda x: np.inf,
    )

    # Denominator
    M0 = density_integral(x, y, mu, std)

    return np.array([[Mxx[0], Mxy[0]], [Myx[0], Myy[0]]]) / M0


com1 = shape_centroid(mu1, std1)
com2 = shape_centroid(mu2, std2)

# Superimpose distributions according to centroids
Q1 = quadrupole(mu1 - com1, std1)
Q2 = quadrupole(mu2 - com2, std2)

print(np.linalg.eig(Q1)[1])
print(np.linalg.eig(Q2)[1])
