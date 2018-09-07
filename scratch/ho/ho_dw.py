import numpy as np
import scipy.linalg
import scipy.sparse
import scipy.sparse.linalg
import matplotlib.pyplot as plt
import time

l = 10
n = 1000
omega = 1
l_x = 2

v_ho = lambda x, omega=1: 0.5 * omega ** 2 * x ** 2
# v_dw = lambda x, omega=1, l_x=1: -omega**2*l_x*np.abs(x) + l_x**2
v = lambda x, omega=1, l_x=1: v_ho(x, omega=omega)  # \
# + v_dw(x, omega=omega, l_x=l_x)

phi_0 = lambda x, omega=1: (omega / np.pi) ** (1 / 4.0) * np.exp(
    -omega * x ** 2 / 2.0
)

x = np.linspace(-l, l, n)
delta_x = x[1] - x[0]


h_diag = 1.0 / (delta_x ** 2) + v(x[1 : n - 1], omega=omega, l_x=l_x)
h_off = -1.0 / float(2 * delta_x ** 2) * np.ones(n - 3)

h = scipy.sparse.diags([h_diag, h_off, h_off], offsets=[0, -1, 1])

t0 = time.time()
epsilon, phi = scipy.sparse.linalg.eigs(h, k=3, which="SM")
t1 = time.time()

print("Time: {0:.3f} sec".format(t1 - t0))

plt.plot(x[1:-1], np.abs(phi[:, 0] / np.sqrt(delta_x)) ** 2)
plt.plot(x[1:-1], np.abs(phi_0(x[1:-1])) ** 2)
plt.show()

__import__("sys").exit()
h = np.zeros((n - 2, n - 2))

for i in range(n - 2):
    h[i, i] = 1.0 / (delta_x ** 2) + v(x[i + 1], omega=omega, l_x=l_x)
    if i + 1 < n - 2:
        h[i + 1, i] = -1.0 / (2 * delta_x ** 2)
        h[i, i + 1] = -1.0 / (2 * delta_x ** 2)


epsilon, phi = scipy.linalg.eigh(h)


plt.plot(x[1:-1], np.abs(phi[:, 0]) ** 2)
plt.show()

plt.plot(x, v(x, omega=omega, l_x=l_x))
plt.show()
