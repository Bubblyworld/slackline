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
