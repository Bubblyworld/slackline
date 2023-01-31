import scipy.optimize as so
import sympy as sp
import numpy as np

def first_order_euler_lagrange(lagrangian):
    """
    Returns a set of 4 coupled first-order ODEs that can be used to integrate
    the solutions to the Euler-Lagrange equations for a given Lagrangian. The
    first two ODEs are for y and n, and the second two are for y_x and n_x.

    The ODEs are returned as a list of numpy functions, where each function
    computes the RHS of d_/dx = f(x, y, n, y_x, n_x, m, g, K).
    """
    # Initialise the symbol variables and functions:
    x = sp.Symbol("x")
    y = sp.Function("y")(x)
    n = sp.Function("n")(x)
    m = sp.Symbol("m")
    g = sp.Symbol("g")
    K = sp.Symbol("K")
    L = lagrangian(x, y, n, m, g, K)

    # The Euler-Lagrange equations can be computed using sympy:
    els = sp.calculus.euler.euler_equations(L, [y, n], [x])
    sol = sp.solve(els, [y.diff(x, 2), n.diff(x, 2)])
    y_xx = sol[y.diff(x, 2)]
    n_xx = sol[n.diff(x, 2)]

    # We introduce auxiliary variables a and b for y_x and n_x respectively:
    a = sp.Function("a")(x)
    b = sp.Function("b")(x)

    # The four ODEs are:
    rhs_y = a
    rhs_n = b
    rhs_a = y_xx.subs({y.diff(x): a, n.diff(x): b})
    rhs_b = n_xx.subs({y.diff(x): a, n.diff(x): b})

    # We lambdify them so that they can be used with numpy arrays:
    fn_y = sp.lambdify([x, y, n, a, b, m, g, K], rhs_y, "numpy")
    fn_n = sp.lambdify([x, y, n, a, b, m, g, K], rhs_n, "numpy")
    fn_a = sp.lambdify([x, y, n, a, b, m, g, K], rhs_a, "numpy")
    fn_b = sp.lambdify([x, y, n, a, b, m, g, K], rhs_b, "numpy")

    return [fn_y, fn_n, fn_a, fn_b]

def mass_boundary_conditions(lagrangian, _x, _m, _g, _K, _M):
    """
    Given the final values of the system at a point mass on its left-hand
    segment of webbing, this function returns the initial conditions of the
    right-hand segment of webbing. The initial conditions are computed using
    the TISE integation trick, by integrating a dirac delta potential for the
    point mass about an infinitesmally small neighbourhood of the mass.

    This kind of trickery is necessary because in general the system has
    discontinuities at point masses, and handling them by simply changing the
    mass of the slackline introduces numerical issues with the integrator.
    """
    x = sp.Symbol("x")
    y = sp.Function("y")(x)
    n = sp.Function("n")(x)
    a_l = sp.Symbol("a_l") # dy/dx on LHS of mass
    b_l = sp.Symbol("b_l") # dn/dx on LHS of mass
    a_r = sp.Symbol("a_r") # dy/dx on RHS of mass
    b_r = sp.Symbol("b_r") # dn/dx the RHS of mass
    m = sp.Symbol("m")
    g = sp.Symbol("g")
    K = sp.Symbol("K")
    M = sp.Symbol("M") # mass of the slackliner

    # First we get dL/da and dL/db:
    L = lagrangian(x, y, n, m, g, K)
    dL_da = sp.diff(L, y.diff(x))
    dL_db = sp.diff(L, n.diff(x))

    # Then we make the substitutions for a_l, b_l, a_r, b_r:
    dL_da_l = dL_da.subs({y.diff(x): a_l, n.diff(x): b_l})
    dL_da_r = dL_da.subs({y.diff(x): a_r, n.diff(x): b_r})
    dL_db_l = dL_db.subs({y.diff(x): a_l, n.diff(x): b_l})
    dL_db_r = dL_db.subs({y.diff(x): a_r, n.diff(x): b_r})

    # We want to solve the following equations for a_r, b_r:
    eq1 = dL_da_r - dL_da_l - M*g
    eq2 = dL_db_r - dL_db_l

    # Substituting for x here might be an issue if y(x), n(x) are around:
    eq1 = eq1.subs({x: _x, m: _m, g: _g, K: _K, M: _M})
    eq2 = eq2.subs({x: _x, m: _m, g: _g, K: _K, M: _M})

    # Recast the equations in a form scipy root solvers expect:
    _fn1 = sp.lambdify([a_l, b_l, a_r, b_r], eq1, "numpy")
    _fn2 = sp.lambdify([a_l, b_l, a_r, b_r], eq2, "numpy")
    fn = lambda a_l, b_l: lambda ab: (
        _fn1(a_l, b_l, ab[0], ab[1]), _fn2(a_l, b_l, ab[0], ab[1])
    )

    # And we're done here:
    return lambda y_x, n_x: _mass_boundary_conditions(fn, y_x, n_x)

def _mass_boundary_conditions(fn, y_x, n_x):
    fn = fn(y_x, n_x)
    root = so.fsolve(fn, [y_x, n_x])
    if np.isclose(fn(root), [0,0], atol=1e-3).all():
        return root
    else:
        return None
