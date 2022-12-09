import numpy as np
import scipy as sp
from .calculus import first_order_euler_lagrange, mass_boundary_conditions
from matplotlib import pyplot as plt


# An out of bounds error for slacklines that are too long:
class SlacklineTooLongError(Exception):
    pass


def integrate(lagrangian, slackline, anchor_tension, anchor_angle,
              masses=[], length_cutoff=10000.0, foel=None):
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
        foel: Optional cache of the first-order Euler-Lagrange equations, for efficiency.
    """
    print(f"Integrating with anchor tension {anchor_tension} and anchor angle {anchor_angle}.")
    if foel is None:
        foel = first_order_euler_lagrange(lagrangian)

    # Remove any masses at the left anchor and warn the user:
    for x, m in masses:
        if x == 0:
            print(f"Warning: removing mass of {m}kg at x=0, as this is the left anchor point.")
    masses = [m for m in masses if m[0] > 0.0]

    # The system of four ODEs we need to integrate:
    fn_y, fn_n, fn_y_x, fn_n_x = foel
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
            L = masses[i][0] - last_x

        # We use a loop because we need to search for the length of the last
        # segment of webbing, which we don't know up front:
        while True:
            if L + last_x > length_cutoff:
                L = length_cutoff - last_x

            _x = np.linspace(last_x, last_x + L, 1000)
            sol = sp.integrate.solve_ivp(
                fn,
                (last_x, last_x + L),
                (y0, n0, y_x0, n_x0),
                method="RK45",
                t_eval=_x,
            )

            # If we've reached the right anchor (y > 0), then we need to trim 
            # the solution and return early:
            if sol.y[0][-1] >= 0.0:
                j = np.argmax(sol.y[0] > 0.0)
                x = np.append(x, sol.t[:j])
                y = np.append(y, sol.y[0][:j])
                n = np.append(n, sol.y[1][:j])
                y_x = np.append(y_x, sol.y[2][:j])
                n_x = np.append(n_x, sol.y[3][:j])

                # Add the right anchor point using a linear interpolation:
                dx = -y[-1] / y_x[-1]
                x = np.append(x, x[-1] + dx)
                y = np.append(y, 0.0)
                n = np.append(n, n[-1] + dx * n_x[-1])
                y_x = np.append(y_x, y_x[-1])
                n_x = np.append(n_x, n_x[-1])

                # If there are outstanding masses, warn the user:
                if i < len(masses):
                    print(f"Warning: {len(masses) - i} masses were ignored as they lie beyond the right anchor point.")

                print(f"Integration complete. Final length: {x[-1]}m.")
                return x, y, n, y_x, n_x

            # If we haven't reached the right anchor point for the last segment,
            # then we need to increase the length of the segment and try again:
            if i == len(masses):
                if L + last_x >= length_cutoff:
                    print(f"Warning: slackline is too long. Stopping at: {L + last_x}m.")
                    raise SlacklineTooLongError(f"Slackline is too long. Maximum length is {length_cutoff}m.")
                L *= 2
                continue

            # Otherwise, we need to apply the boundary condition at the mass
            # and continue integrating the next segment:
            x = np.append(x, sol.t)
            y = np.append(y, sol.y[0])
            n = np.append(n, sol.y[1])
            y_x = np.append(y_x, sol.y[2])
            n_x = np.append(n_x, sol.y[3])
            last_x = x[-1]
            print(f"last_y: {y[-1]}")
            y0, n0, y_x0, n_x0 = mass_boundary_conditions(
                lagrangian, last_x, y[-1], n[-1], y_x[-1], n_x[-1],
                slackline.m, slackline.g, slackline.K, masses[i][1],
            )
            break

    print(f"Integration complete. Final length: {x[-1]}m.")
    return x, y, n, y_x, n_x


def integrate_length_tension(lagrangian, slackline, gap_length, anchor_tension,
                             masses=[], foel=None):
    """
    Integrator that solves for a slackline of a given gap length and standing
    anchor tension by binary searching the core integrator.
    """
    print(f"Integrating for gap length: {gap_length}m, anchor tension: {anchor_tension}N.")
    if foel is None:
        foel = first_order_euler_lagrange(lagrangian)

    a = 0.001 # lower bound
    b = np.pi / 4 # upper bound
    while True:
        anchor_angle = (a + b) / 2
        try:
            x, y, n, y_x, n_x = integrate(
                lagrangian, slackline, anchor_tension, anchor_angle,
                masses=masses,
                length_cutoff=gap_length,
                foel=foel,
            )
        except SlacklineTooLongError:
            b = anchor_angle
            continue
        if np.abs(x[-1] - gap_length) < 1e-1: # don't have to be too precise
            return x, y, n, y_x, n_x
        elif x[-1] > gap_length:
            b = anchor_angle
        else:
            a = anchor_angle

def integrate_natural_length(lagrangian, slackline, gap_length, natural_length,
                             masses=[], foel=None):
    """
    Integrator that solves for a slackline of a given gap length and natural
    (i.e. untensioned) slackline length. These parameters can be obtained from
    the length_tension integrator, and then used to solve for tension once the
    slackliners have been positioned.
    """
    print(f"Integrating for natural length: {natural_length}m, gap length: {gap_length}m.")
    if foel is None:
        foel = first_order_euler_lagrange(lagrangian)

    # We can binary search on the anchor tension:
    a = natural_length * slackline.m * slackline.g # lower bound
    b = 50000.0 # upper bound
    while True:
        anchor_tension = (a + b) / 2
        x, y, n, y_x, n_x = integrate_length_tension(
            lagrangian, slackline, gap_length, anchor_tension,
            masses=masses,
            foel=foel,
        )
        if np.abs(n[-1] - natural_length) < 1e-1: # don't have to be too precise
            return x, y, n, y_x, n_x, anchor_tension
        elif n[-1] > natural_length:
            a = anchor_tension
        else:
            b = anchor_tension
