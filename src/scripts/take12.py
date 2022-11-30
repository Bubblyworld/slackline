import matplotlib.pyplot as plt
import numpy as np
import sympy as sp

# Use unicode to print sympy equations:
sp.init_printing(use_unicode=True)

# Constants:
dx = 0.1 # m (step size)
g = 9.81 # m / s^2 (gravitational acceleration)
K = 2500*100.0 # N (newtons per 100% stretch in rope)
m = 0.044 # kg / m (mass density of rope in natural units)

# Initial conditions:
L = 100.0 # m (gap length)
T0 = 1000.0 # N (initial tension in rope)
theta0 = 0.03 # rads (initial angle of webbing below horizontal)

def model_diagnostics():
    """
    Uses SymPy to check a bunch of the algebra behind the model.
    """
    print("Starting model diagnostics.")

    # The two functions we are trying to find are y(x) and n(x), which are the
    # height of the webbing at a point x and the natural length of the webbing
    # at a point x, respectively.
    x = sp.Symbol("x")
    print("   ...creating independent variable x.")
    y = sp.Function("y")(x)
    n = sp.Function("n")(x)
    print("   ...creating dependent functions y(x) and n(x).")

    # The model constants are the following:
    #   g: gravitational acceleration
    #   K: newtons per 100% stretch in rope
    #   m: mass density of rope in units of n
    g = sp.Symbol("g")
    K = sp.Symbol("K")
    m = sp.Symbol("m")

    # The potential energy due to gravity for a segment of webbing is:
    #   dU_g = m*g*y dn = m*g*y*n_x dx
    # ...where we denote the derivative of n with respect to x by n_x.
    U_g = m*g*y*n.diff(x)
    print()
    print("   ...potential energy due to gravity is given by:")
    print(U_g)

    # The potential energy due to the tension in the rope is:
    #   dU_T = [ K/2*(2 - n_x)*(1 + y_x^2) - K/2*y_x^2 + K/2*n_x ] dx
    # ...where we denote the derivative of y with respect to x by y_x.
    U_T = K/2*(1/n.diff(x))*(1 + y.diff(x)**2) - K/2*y.diff(x)**2 + K/2*n.diff(x)
    print()
    print("   ...potential energy due to tension in the rope is given by:")
    print(U_T)

    # The Lagrangian is the sum of the potential energies:
    L = U_g + U_T
    print()
    print("   ...the Lagrangian is given by:")
    print(L)

    # The Euler-Lagrange equations are determined by the derivatives of the
    # Lagrangian with respect to y, y_x, n and n_x:
    dL_dy = sp.diff(L, y)
    dL_dy_x = sp.diff(L, y.diff(x))
    dL_dn = sp.diff(L, n)
    dL_dn_x = sp.diff(L, n.diff(x))
    print()
    print("   ...the components of the Euler-Lagrange equations are:")
    print("dL/dy   = {}".format(dL_dy))
    print("dL/dy_x = {}".format(dL_dy_x))
    print("dL/dn   = {}".format(dL_dn))
    print("dL/dn_x = {}".format(dL_dn_x))

    # We can solve the Euler-Lagrange equations directly in sympy:
    print()
    print("   ...solving the Euler-Lagrange equations directly in sympy:")
    els = sp.calculus.euler.euler_equations(L, [y, n], [x])
    print("[0]: {}".format(els[0]))
    print("[1]: {}".format(els[1]))

    # Use sympy to solve the resulting equations:
    print()
    print("   ...solving the Euler-Lagrange equations using sympy:")
    sol = sp.solve(els, [y.diff(x, 2), n.diff(x, 2)])
    y_xx = sol[y.diff(x, 2)]
    n_xx = sol[n.diff(x, 2)]
    print("y_xx = {}".format(y_xx))
    print("n_xx = {}".format(n_xx))

    # Initial conditions y(0), n(0), y_x(0), n_x(0) are:
    y0 = 0.0
    n0 = 0.0
    y_x0 = -0.001
    n_x0 = 0.999

    # Find initial y_xx(0) and n_xx(0) by substituting initial conditions:
    y_xx0 = y_xx.subs([(y.diff(x), y_x0), (n.diff(x), n_x0), (y, y0), (n, n0)])
    n_xx0 = n_xx.subs([(y.diff(x), y_x0), (n.diff(x), n_x0), (y, y0), (n, n0)])
    print()
    print("   ...initial conditions are:")
    print("y(0)   = {}".format(y0))
    print("n(0)   = {}".format(n0))
    print("y_x(0) = {}".format(y_x0))
    print("n_x(0) = {}".format(n_x0))
    print()
    print("   ...initial values of y_xx(0) and n_xx(0) are:")
    print("y_xx(0) = {}".format(y_xx0))
    print("n_xx(0) = {}".format(n_xx0))

    # Subsitute known values for K, g, m:
    y_xx0 = y_xx0.subs([(K, 2500*100.0), (g, 9.81), (m, 0.044)])
    n_xx0 = n_xx0.subs([(K, 2500*100.0), (g, 9.81), (m, 0.044)])
    print()
    print("   ...substituting known values for K, g, m:")
    print("y_xx(0) = {}".format(y_xx0))
    print("n_xx(0) = {}".format(n_xx0))

def main():
    model_diagnostics()
    return

    # Variables:
    x = [0.0]
    y = [0.0]
    n = [0.0]

    # Compute initial conditions y_x0, n_x0 of the model using the simplified
    # initial conditions T0, theta0:
    y_x = [np.tan(-theta0)]
    n_x = [np.sqrt(1 + y_x[0]**2) / (T0/K + 1)]

    # Compute the constants of motion:
    A_1 = -K/2 * y_x[0]**2
    A_2 = K * y_x[0] * (1 - n_x[0])

    # Integrate the model, which is given by the following:
    i = 0
    while y_x[-1] < -1e-4 and n_x[-1] < 1.1: # the model has problems switching direction
        x.append(x[i] + dx)
        y.append(y[i] + dx * y_x[i])
        n.append(n[i] + dx * n_x[i])
        y_x.append(-np.sqrt((m * g * y[i+1] - A_1) * 2 / K))
        n_x.append(1 - (A_2 + m * g * n[i+1]) / K / y_x[i+1])
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
    axes[0, 0].set_title("y")
    axes[0, 1].plot(x, y_x)
    axes[0, 1].set_title("y_x")
    axes[1, 0].plot(x, n)
    axes[1, 0].set_title("n")
    axes[1, 1].plot(x, n_x)
    axes[1, 1].set_title("n_x")
    plt.show()
