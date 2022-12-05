import sympy as sp

def ideal(x, y, n, m, g, K):
    """
    A Lagrangian for an ideal elastic webbing with the following parameters:
      x: horizontal distance along the gap
      y(x): vertical drop of webbing
      n(x): unstretched length of the webbing up to point x
      m: mass per meter of the webbing
      g: gravitational acceleration
      K: newtons per 100% stretch in rope
    """
    y_x = y.diff(x)
    n_x = n.diff(x)
    gravity = m*g*y*n_x
    tension = K/2*(1 + y_x**2)/n_x - K*sp.sqrt(1 + y_x**2) + K/2*n_x
    return gravity + tension

def first_approximation(x, y, n, m, g, K):
    """
    An apprimation of the ideal Lagrangian that makes the following assumptions:
        1. sqrt(1 + y_x^2) = 1 + y_x^2/2
    """
    return ideal(x, y, n, m, g, K).subs({
        sp.sqrt(1 + y.diff(x)**2): 1 + y.diff(x)**2/2,
    })
