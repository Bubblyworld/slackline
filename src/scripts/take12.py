import matplotlib.pyplot as plt
import numpy as np
from findiff import FinDiff

# Constants:
dx = 0.001 # m (step size)
g = 9.81 # m / s^2 (gravitational acceleration)
K = 2500*100.0 # N (newtons per 100% stretch in rope)
m = 0.044 # kg / m (mass density of rope in natural units)

# Initial conditions:
L = 100.0 # m (gap length)
T0 = 1000.0 # N (initial tension in rope)
theta0 = 0.03 # rads (initial angle of webbing below horizontal)

def main():
    # Variables:
    x = [0.0]
    y = [0.0]
    n = [0.0]

    # Compute initial conditions y_x0, n_x0 of the model using the simplified
    # initial conditions T0, theta0:
    y_x = [np.tan(-theta0)]
    n_x = [np.sqrt(1 + y_x[0]**2) / (T0/K + 1)]

    # Integrate the model, which is given by the following:
    i = 0
    while y_x[-1] < 0.0: # the model has problems switching direction
        x.append(x[i] + dx)
        y.append(y[i] + dx * y_x[i])
        n.append(n[i] + dx * n_x[i])
        y_xx = m * g * y[i] / K / y_x[i]
        y_x.append(y_x[i] + dx * y_xx)
        n_xx = (K * y_xx * (1 - n_x[i]) - m * g * n_x[i]) / K / y_x[i]
        n_x.append(n_x[i] + dx * n_xx)
        i += 1

    # Convert everything to numpy and compute some utility variables:
    x = np.array(x)
    y = np.array(y)
    n = np.array(n)
    y_x = np.array(y_x)
    n_x = np.array(n_x)
    T = K * (np.sqrt(1 + y_x**2) / n_x - 1) # something weird here, negative?
    l = np.cumsum(np.sqrt(1 + y_x**2) * dx)

    # Plot the results as a 2x2 grid of plots: y, n, l, T
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    axes[0, 0].plot(x, y)
    # set y range to [-100,0]
    axes[0, 0].set_ylim(-100, 0)
    axes[0, 0].set_title("y")
    axes[0, 1].plot(x, y_x)
    axes[0, 1].set_title("y_x")
    axes[1, 0].plot(x, l)
    axes[1, 0].set_title("l")
    axes[1, 1].plot(x, T)
    axes[1, 1].set_title("T")
    plt.show()
