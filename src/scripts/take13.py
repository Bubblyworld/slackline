import matplotlib.pyplot as plt
import numpy as np
import sympy as sp
from findiff import FinDiff

# Setup:
sp.init_printing()

# Constants:
dx = 0.1 # m (step size)
g = 9.81 # m / s^2 (gravitational acceleration)
K = 2500*100.0 # N (newtons per 100% stretch in rope)
m = 0.044 # kg / m (mass density of rope in natural units)

# Initial conditions:
X = 100.0 # m (gap length)
T0 = 1000.0 # N (initial tension in rope)
theta0 = 0.03 # rads (initial angle of webbing below horizontal)

################################################################################
# The two functions we are trying to find are y(x) and n(x), which are the     #
# height of the webbing at a point x and the natural length of the webbing     #
# at a point x, respectively. To solve for y(x) and n(x), we apply variational #
# calculus to an appropriate Lagrangian, something sympy can do for us.        #
################################################################################

def lagrangian(x, y, n, g, m, K):
    """
    A Lagrangian with no approximations.
    """
    y_x = y.diff(x)
    n_x = n.diff(x)
    gravity = m*g*y*n.diff(x)
    tension = K/2*(1 + y_x**2)/n_x - K*sp.sqrt(1 + y_x**2) + K/2*n_x
    return gravity + tension

def lagrangian_approx_1(x, y, n, g, m, K):
    """
    A Lagrangian that uses the following second-order approximation:
        sqrt(1 + y_x^2) ~ 1 + y_x^2/2
    ...which is only valid in the regime y_x ~ 0.
    """
    return lagrangian(x, y, n, g, m, K).subs({
        sp.sqrt(1 + y.diff(x)**2): 1 + y.diff(x)**2/2,
    })

################################################################################
# To integrate y(x) and n(x) numerically, we use the Euler-Lagrange equations. #
# derived from the Lagrangian. These are second-order ODEs that give y_xx and  #
# n_xx as functions of n, y, n_x and y_x.                                      #
################################################################################

def integrate(lagrangian):
    # Initialise the symbol variables and functions:
    _x = sp.Symbol("x")
    _y = sp.Function("y")(_x)
    _n = sp.Function("n")(_x)
    _m = sp.Symbol("m")
    _g = sp.Symbol("g")
    _K = sp.Symbol("K")
    _L = lagrangian(_x, _y, _n, _g, _m, _K)
    sp.pprint(_L)

    # The Euler-Lagrange equations can be computed using sympy:
    els = sp.calculus.euler.euler_equations(_L, [_y, _n], [_x])
    sol = sp.solve(els, [_y.diff(_x, 2), _n.diff(_x, 2)])
    _y_xx = sol[_y.diff(_x, 2)]
    _n_xx = sol[_n.diff(_x, 2)]

    # We can substitute known values for m, g, K at this point:
    _y_xx = _y_xx.subs({_m: m, _g: g, _K: K})
    _n_xx = _n_xx.subs({_m: m, _g: g, _K: K})

    # And finally, we simplify:
    _y_xx = sp.simplify(_y_xx)
    _n_xx = sp.simplify(_n_xx)

    # Print them out:
    print("EULER-LAGRANGE EQUATIONS:")
    print("y_xx =")
    sp.pprint(_y_xx)
    print()
    print("n_xx =")
    sp.pprint(_n_xx)
    print()

    # Initialise the numerical variables and functions using numpy arrays:
    x = np.array([0.0])
    y = np.array([0.0])
    n = np.array([0.0])

    # Compute initial conditions y_x0, n_x0 of the model using the simplified
    # initial conditions T0, theta0:
    y_x = np.array([np.tan(-theta0)])
    n_x = np.array([np.sqrt(1 + y_x[0]**2) / (T0/K + 1)])

    # Print them out:
    print("INITIAL CONDITIONS:")
    print("y_x0 = {}".format(y_x[0]))
    print("n_x0 = {}".format(n_x[0]))
    print()

    # Actually integrate the model:
    print("INTEGRATING UNTIL y(x) = 0 AGAIN:")
    i=0
    while y[-1] <= 0: # don't run forever
        # We can get y_xx and n_xx from the Euler-Lagrange equations:
        subs = {
            _y.diff(_x): y_x[i],
            _n.diff(_x): n_x[i],
            _y: y[i],
            _n: n[i],
            _x: x[i],
        }
        y_xx = _y_xx.subs(subs)
        n_xx = _n_xx.subs(subs)
        y_xx = sp.lambdify((), y_xx)()
        n_xx = sp.lambdify((), n_xx)()

        # We can then integrate everything directly:
        y_x = np.append(y_x, y_x[i] + y_xx*dx)
        n_x = np.append(n_x, n_x[i] + n_xx*dx)
        y = np.append(y, y[i] + y_x[i]*dx)
        n = np.append(n, n[i] + n_x[i]*dx)
        x = np.append(x, x[i] + dx)
        i += 1

        # Spit something out so we know it's still working:
        if i % 1000 == 0:
            print("step {:5d}: x = {:10.5f}, y = {:10.5f}, n = {:10.5f}".format(
                i, x[i], y[i], n[i]))
    print("Integration complete after {} steps.".format(i))
    print("DONE.")

    # Compute arclength and tension:
    l = np.cumsum(np.sqrt(1 + y_x**2))*dx
    T = K*(np.sqrt(1 + y_x**2) / n_x - 1)

    # Return the results:
    return x, y, n, y_x, n_x, l, T

def main():
    x1, y1, n1, y_x1, n_x1, l1, T1 = integrate(lagrangian)
    x2, y2, n2, y_x2, n_x2, l2, T2 = integrate(lagrangian_approx_1)

    # Plot y(x) for all Lagrangians:
    plt.figure()
    plt.plot(x1, y1, label="Lagrangian")
    plt.plot(x2, y2, label="Lagrangian (approx 1)")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.legend()
    plt.show()
