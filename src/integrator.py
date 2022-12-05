import numpy as np
import scipy as sp
from .calculus import first_order_euler_lagrange
from matplotlib import pyplot as plt


# An out of bounds error for slacklines that are too long:
class SlacklineTooLongError(Exception):
    pass


def integrate(lagrangian, slackline, anchor_tension, anchor_angle, length_cutoff=10000.0):
    """
    Integrates the curve of a rigged slackline with given anchor tension, and
    angle. This can be binary searched to do more interesting things, like
    finding the curve for a slackline of a given gap length.
    """
    # The system of four ODEs we need to integrate:
    fn_y, fn_n, fn_y_x, fn_n_x = first_order_euler_lagrange(lagrangian)

    # Initial conditions for each variable:
    y0 = 0.0
    n0 = 0.0
    y_x0 = np.tan(-anchor_angle)
    n_x0 = np.sqrt(1 + y_x0**2) / (anchor_tension / slackline.K + 1)

    # We're not sure how long the slackline will be, so we need to search for
    # that, starting with a guess of X=1000m.
    L = 1000.0
    while L < length_cutoff:
        x = np.linspace(0, L, 2000)
        def fn(x, y):
            return np.array([
                fn_y(x, y[0], y[1], y[2], y[3], slackline.m, slackline.g, slackline.K),
                fn_n(x, y[0], y[1], y[2], y[3], slackline.m, slackline.g, slackline.K),
                fn_y_x(x, y[0], y[1], y[2], y[3], slackline.m, slackline.g, slackline.K),
                fn_n_x(x, y[0], y[1], y[2], y[3], slackline.m, slackline.g, slackline.K),
            ])
        sol = sp.integrate.solve_ivp(
            fn,
            (0, L),
            (y0, n0, y_x0, n_x0),
            method="BDF",
            t_eval=x,
        )

        # If the curve hasn't reached the other anchor, we need to increase
        # the length of the curve we're integrating by double:
        if sol.y[0][-1] < 0:
            L *= 2
            continue

        # If the curve has reached the other anchor, we need to cap the arrays
        # to the correct length:
        x = sol.t
        y = sol.y[0]
        n = sol.y[1]
        y_x = sol.y[2]
        n_x = sol.y[3]
        i = np.argmax(y > 0)
        return x[:i], y[:i], n[:i], y_x[:i], n_x[:i]

    # If we reach this point, we've failed to find the end of the slackline:
    raise SlacklineTooLongError("The slackline is too long to integrate.")
