import sympy as sp

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

def mass_boundary_conditions(lagrangian, x, y, n, y_x, n_x, m, g, K, mass):
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
    # We assume continuity of x, y(x) and n(x), so the initial values for
    # these variables is simply the final values of the previous segment:
    y0 = y
    n0 = n

    # To compute the initial values of y_x and n_x, we can solve the following
    # two equations obtained by integrating a delta potential over [x-e, x+e]:
    #    dL/dy_x(right) - dL/dy_x(left) = m*g
    #    dL/dn_x(right) - dL/dn_x(left) = 0
    # where L is the Lagrangian of the system.
    _x = sp.Symbol("x")
    _y = sp.Function("y")(_x)
    _n = sp.Function("n")(_x)
    _m = sp.Symbol("m")
    _g = sp.Symbol("g")
    _K = sp.Symbol("K")
    _y_x = _y.diff(_x)
    _n_x = _n.diff(_x)
    L = lagrangian(_x, _y, _n, _m, _g, _K)
    dL_dy_x = sp.diff(L, _y_x)
    dL_dn_x = sp.diff(L, _n_x)

    # Sympy has issues with substitutions, so we introduce new symbols for y_x, n_x:
    a = sp.Symbol("a")
    b = sp.Symbol("b")
    dL_dy_x = dL_dy_x.subs({_y_x: a, _n_x: b})
    dL_dn_x = dL_dn_x.subs({_y_x: a, _n_x: b})

    # For the left-side values, we substitute the known variables:
    dL_dy_x_left = dL_dy_x.subs({_x: x, _y: y, _n: n, _m: m, _g: g, _K: K, a: y_x, b: n_x})
    dL_dn_x_left = dL_dn_x.subs({_x: x, _y: y, _n: n, _m: m, _g: g, _K: K, a: y_x, b: n_x})

    # For the right-side values, we substitute everything but a, b:
    dL_dy_x_right = dL_dy_x.subs({_x: x, _y: y, _n: n, _m: m, _g: g, _K: K})
    dL_dn_x_right = dL_dn_x.subs({_x: x, _y: y, _n: n, _m: m, _g: g, _K: K})

    # Construct and solve the two equations:
    eq1 = sp.Eq(dL_dy_x_right - dL_dy_x_left, mass*g)
    eq2 = sp.Eq(dL_dn_x_right - dL_dn_x_left, 0)
    sol = sp.solve([dL_dy_x_right - dL_dy_x_left - mass*g, dL_dn_x_right - dL_dn_x_left], [a, b])

    # The true solutions must have a, b both real with b>0:
    for (a, b) in sol:
        if a.is_real and b.is_real and b > 0:
            return [y0, n0, a, b]

    # If no solution is found, raise an exception:
    raise Exception("No solution found for initial conditions!")
