import partitionunity as pu
import numpy as np
import time
import matplotlib.pyplot as plt

# Define the test function
def franke(X):
    return 0.75 * np.exp(-(9*X[:, 0]-2)**2/4 - (9*X[:, 1]-2)**2/4) + \
           0.75 * np.exp(-(9*X[:, 0]+1)**2/49 - (9*X[:, 1]+1) / 10) + 0.5 * \
           np.exp(-(9*X[:, 0]-7)**2/4 - (9*X[:, 1]-3)**2/4)-0.2 * \
           np.exp(-(9*X[:, 0]-4)**2 - (9*X[:, 1]-7)**2)


# Define the space dimension and the number of interpolation data
d = 2
p = np.append(2, np.arange(2, 7))
n = np.floor(np.sqrt(7**p)).astype(int)


# Define the kernel and parameter
def Phi(eps, r):
    return (1 + eps * r) * np.exp(-eps * r)


epsilon = 1

# Define the weights for PU
def weight(e, r):
    return np.multiply(np.power(np.fmax(1-(e*r), 0*(e*r)), 4), (4*(e*r)+1))


# Define s_d^d equally spaced test data
s_d = 60
bar_x = np.array(np.meshgrid(np.linspace(0, 1, s_d), np.linspace(0, 1, s_d))).T.reshape(-1, d)
t = []
MAE = []
for i in range(len(p)):

    # Define the interpolation data
    x = np.array(np.meshgrid(np.linspace(0, 1, n[i]), np.linspace(0, 1, n[i]))).T.reshape(-1, d)

    # Define the function values
    y = franke(x)

    m_d = np.floor((n[i] ** 2) ** (1 / d) / 2)
    tm = time.time()

    # Compute the PU interpolant
    Pf = pu.PU(x, y, bar_x, m_d, weight, Phi, epsilon)

    t.append(time.time() - tm)
    MAE.append(np.max(abs(franke(bar_x) - Pf)))

# Display the results
plt.loglog(n[1:]**2, MAE[1:], 'o-', )
plt.title("MAE vs number of points")
plt.xlabel("n")
plt.ylabel("MAE")
plt.grid(True, which="both", linestyle='--')
plt.show()

plt.loglog(n[1:]**2, t[1:], 'o-')
plt.title("CPU time vs number of points")
plt.xlabel("n")
plt.ylabel("CPU(s)")
plt.grid(True, which="both", linestyle='--')
plt.show()

