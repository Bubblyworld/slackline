import sympy as sp
import numpy as np

###############################################################################
# Lagrangians are defined in terms of x, y(x), n(x), m(n), g, K.               #
###############################################################################

def lagrangian(x, y, n, m, g, K):
    """
    A Lagrangian for a realistic elastic webbing with no approximations.
    """
    y_x = y.diff(x)
    n_x = n.diff(x)
    gravity = m*g*y*n.diff(x)
    tension = K/2*(1 + y_x**2)/n_x - K*sp.sqrt(1 + y_x**2) + K/2*n_x
    return gravity + tension

def lagrangian_approx_1(x, y, n, m, g, K):
    """
    A Lagrangian for a realistic elastic webbing with the following
    second-order approximation:
        sqrt(1 + y_x^2) ~ 1 + y_x^2/2
    ...which is only valid in the regime y_x ~ 0.
    """
    return lagrangian(x, y, n, m, g, K).subs({
        sp.sqrt(1 + y.diff(x)**2): 1 + y.diff(x)**2/2,
    })

def lagrangian_dd(x, y, n, m, g, K):
    """
    A Lagrangian that includes a dirac-delta term at n=0 for a mass of weight
    W. This is used to simulate the mass of the webbing itself.
    """
    return lagrangian_approx_1(x, y, n, m, g, K) + m*sp.DiracDelta(n)*y

###############################################################################
# Run the script.                                                             #
###############################################################################

def main():
    """
    The idea with this script is to derive the equations of motion for various
    kinds of webbing lagrangians. We can then convert these equations into a
    system of 1st order ODEs and integrate them with scypy, hopefully with a
    fast, accurate numerical integrator.
    """
    # Create symbolic variables:
    g = sp.Symbol('g')
    K = sp.Symbol('K')
    x = sp.Symbol('x')
    y = sp.Function('y')(x)
    n = sp.Function('n')(x)
    m = sp.Function('m')(n)
    L = lagrangian_dd(x, y, n, m, g, K)

    # Pretty print the Lagrangian:
    print("LAGRANGIAN:")
    sp.pprint(L)
    print()

    # Compute the Euler-Lagrange equations:
    els = sp.calculus.euler.euler_equations(L, [y, n], [x])
    sol = sp.solve(els, [y.diff(x, 2), n.diff(x, 2)])

    # Print the solutions:
    print("EULER-LAGRANGE EQUATIONS:")
    print("y_xx =")
    sp.pprint(sol[y.diff(x, 2)])
    print()
    print("n_xx =")
    sp.pprint(sol[n.diff(x, 2)])
    print()

    # Convert the two 2nd order ODEs into a system of 1st order ODEs by
    # introducing auxilliary variables a(x) = y.diff(x) and b(x) = n.diff(x):
    a = sp.Function('a')(x)
    b = sp.Function('b')(x)
    y_xx = sol[y.diff(x, 2)]
    n_xx = sol[n.diff(x, 2)]
    eq_a = sp.Eq(a.diff(x), y_xx.subs({y.diff(x): a, n.diff(x): b}))
    eq_b = sp.Eq(b.diff(x), n_xx.subs({y.diff(x): a, n.diff(x): b}))
    eq_y = sp.Eq(y.diff(x), a)
    eq_n = sp.Eq(n.diff(x), b)
    print("SYSTEM OF 1ST ORDER ODEs:")
    sp.pprint(eq_a)
    sp.pprint(eq_b)
    sp.pprint(eq_y)
    sp.pprint(eq_n)
    print()
