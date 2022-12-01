import matplotlib.pyplot as plt
import numpy as np
import sympy as sp
import scipy
from findiff import FinDiff

# Setup:
sp.init_printing()

# Constants:
dx = 0.1 # m (step size)
g = 9.81 # m / s^2 (gravitational acceleration)
K = 2500*100.0 # N (newtons per 100% stretch in rope)
m = 0.044 # kg / m (mass density of rope in natural units)

# Constants for the boundary-value problem:
X = 100.0 # m (gap length)
N = 99.0 # m (natural length)

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
# n_xx as functions of n, y, n_x and y_x. We can then convert these to first-  #
# order ODEs by defining a = dy/dx and b = dn/dx.                              #
################################################################################

def integrate(lagrangian):
    # Initialise the symbol variables and functions:
    _x = sp.Symbol("x")
    _y = sp.Function("y")(_x)
    _n = sp.Function("n")(_x)
    _m = sp.Function("m")(_n)
    _g = sp.Symbol("g")
    _K = sp.Symbol("K")
    _L = lagrangian(_x, _y, _n, _g, _m, _K)

    # Spit out the Lagrangian so we know what we're working with:
    print("LAGRANGIAN:")
    sp.pprint(_L)
    print()

    # The Euler-Lagrange equations can be computed using sympy:
    els = sp.calculus.euler.euler_equations(_L, [_y, _n], [_x])
    sol = sp.solve(els, [_y.diff(_x, 2), _n.diff(_x, 2)])
    _y_xx = sol[_y.diff(_x, 2)]
    _n_xx = sol[_n.diff(_x, 2)]

    # Now we can introduce the auxilliary variables a = dy/dx and b = dn/dx:
    _a = sp.Function("a")(_x)
    _b = sp.Function("b")(_x)
    eq_a = _y_xx.subs({ _y.diff(_x): _a, _n.diff(_x): _b })
    eq_b = _n_xx.subs({ _y.diff(_x): _a, _n.diff(_x): _b })
    eq_y = _a
    eq_n = _b

    # Spit out the system of first-order ODEs:
    print("SYSTEM OF FIRST-ORDER ODEs:")
    print("dy/dx =")
    sp.pprint(eq_y)
    print()
    print("dn/dx =")
    sp.pprint(eq_n)
    print()
    print("da/dx =")
    sp.pprint(eq_a)
    print()
    print("db/dx =")
    sp.pprint(eq_b)
    print()

    # Lambdify the equations so we can use them in a numpy integrator:
    fn_y = sp.lambdify((_x, _y, _n, _a, _b, _m, _g, _K), eq_y, "numpy")
    fn_n = sp.lambdify((_x, _y, _n, _a, _b, _m, _g, _K), eq_n, "numpy")
    fn_a = sp.lambdify((_x, _y, _n, _a, _b, _m, _g, _K), eq_a, "numpy")
    fn_b = sp.lambdify((_x, _y, _n, _a, _b, _m, _g, _K), eq_b, "numpy")

    # Our boundary conditions are given by y(0) = n(0) = 0, y(X) = 0, and
    # n(X) = N. We can use these to define the initial conditions:
    def bc(ya, yb):
        return np.array([ya[0], ya[1], yb[0], yb[1] - N])

    # Our initial guesses for the solution are y(x) = 0, n(x) = N*x/X,
    # a(x) = 0, and b(x) = N/X:
    x = np.arange(0, X, 0.1)
    guesses = np.zeros((4, x.size))
    guesses[1] = N*x/X
    guesses[3] = N/X * np.ones(x.size)

    # Use scipy.integrate.solve_bvp to solve this boundary-value problem:
    print("SOLVING...")
    def fn(x, y):
        return np.array([
            fn_y(x, y[0], y[1], y[2], y[3], m, g, K),
            fn_n(x, y[0], y[1], y[2], y[3], m, g, K),
            fn_a(x, y[0], y[1], y[2], y[3], m, g, K),
            fn_b(x, y[0], y[1], y[2], y[3], m, g, K),
        ])
    sol = scipy.integrate.solve_bvp(fn, bc, x, guesses)
    print("...and we're done.")
    print()

    # Unpack the integrated curves:
    y = sol.y[0]
    n = sol.y[1]
    y_x = sol.y[2]
    n_x = sol.y[3]

    # Return the solutions with some auxilliary variables:
    l = np.cumsum(np.sqrt(1 + y_x**2))*dx
    T = K*(np.sqrt(1 + y_x**2) / n_x - 1)
    return x, y, n, y_x, n_x, l, T

def main():
    x, y, n, y_x, n_x, l, T = integrate(lagrangian)

    # Plot y against x for each integration with an appropriate label:
    plt.plot(x, y)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.show()

