import numpy as np
import scipy as sp
from .calculus import first_order_euler_lagrange
from matplotlib import pyplot as plt


# An out of bounds error for slacklines that are too long:
class SlacklineTooLongError(Exception):
    pass


def integrate(lagrangian, slackline, anchor_tension, anchor_angle, 
              masses=[], length_cutoff=10000.0):
    """
    Integrates the curve of a rigged slackline with given anchor tension, and
    angle. This can be binary searched to do more interesting things, like
    finding the curve for a slackline of a given gap length.

    Inputs:
        lagrangian: A sympy expression for the Lagrangian of the system.
        slackline: A Slackline object that describes the webbing.
        anchor_tension: The tension in the webbing at the left anchor point.
        anchor_angle: The angle below horizontal of the webbing at the left anchor point.
        masses: A list of (x, m) tuples, where x is the position of the mass, and m is the mass. Defaults to the empty list.
        length_cutoff: The maximum length of the slackline. Defaults to 10000m.
    """
    # The system of four ODEs we need to integrate:
    fn_y, fn_n, fn_y_x, fn_n_x = first_order_euler_lagrange(lagrangian)
    def fn(x, y):
        return np.array([
            fn_y(x, y[0], y[1], y[2], y[3], slackline.m, slackline.g, slackline.K),
            fn_n(x, y[0], y[1], y[2], y[3], slackline.m, slackline.g, slackline.K),
            fn_y_x(x, y[0], y[1], y[2], y[3], slackline.m, slackline.g, slackline.K),
            fn_n_x(x, y[0], y[1], y[2], y[3], slackline.m, slackline.g, slackline.K),
        ])

    # Initial conditions for each variable at the left anchor:
    y0 = 0.0
    n0 = 0.0
    y_x0 = np.tan(-anchor_angle)
    n_x0 = np.sqrt(1 + y_x0**2) / (anchor_tension / slackline.K + 1)

    # We treat masses as point particles, which introduce discontinuities in 
    # the system and require us to integrate the segments of webbing between
    # them independently.
    x = np.array([])
    y = np.array([])
    n = np.array([])
    y_x = np.array([])
    n_x = np.array([])
    last_x = 0.0
    for i in range(len(masses) + 1):
        # If this is the last segment, we integrate until we reach y > 0,
        # which represents the right anchor point. Otherwise, we integrate
        # until we reach the next mass:
        if i == len(masses):
            L = 1000.0
        else:
            L = masses[i][0]

        # We use a loop because we need to search for the length of the last
        # segment of webbing, which we don't know up front:
        while True:
            if L + last_x > length_cutoff:
                raise SlacklineTooLongError("The slackline is too long to integrate.")

            _x = np.linspace(last_x, last_x + L, 1000)
            sol = sp.integrate.solve_ivp(
                fn,
                (last_x, last_x + L),
                (y0, n0, y_x0, n_x0),
                method="BDF",
                t_eval=_x,
            )

            # Path for the last segment:
            if i == len(masses):
                # If we haven't reached the right anchor point, we need to
                # increase the length of the segment:
                if sol.y[0][-1] < 0.0:
                    L *= 2
                    continue

                # Otherwise we need to trim the solution to the right anchor:
                j = np.argmax(sol.y[0] > 0.0)
                x = np.append(x, sol.t[:j])
                y = np.append(y, sol.y[0][:j])
                n = np.append(n, sol.y[1][:j])
                y_x = np.append(y_x, sol.y[2][:j])
                n_x = np.append(n_x, sol.y[3][:j])
                break

            # Path for an initial segment:
            x = np.append(x, sol.t)
            y = np.append(y, sol.y[0])
            n = np.append(n, sol.y[1])
            y_x = np.append(y_x, sol.y[2])
            n_x = np.append(n_x, sol.y[3])
            last_x = x[-1]

            # TODO: update the initial conditions based on dirac delta stuff
            y0 = y[-1]
            n0 = n[-1]
            y_x0 = y_x[-1]
            n_x0 = n_x[-1]
            break

    # ...and we're done here:
    return x, y, n, y_x, n_x
