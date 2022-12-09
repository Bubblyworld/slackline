import matplotlib.pyplot as plt
import numpy as np
from .calculus import first_order_euler_lagrange, mass_boundary_conditions
from .lagrangians import ideal, first_approximation
from .slacklines import dyneemite_pro
from .integrator import integrate, integrate_length_tension, integrate_natural_length

def main():
    x, y, n, y_x, n_x = integrate_length_tension(
        ideal, dyneemite_pro, 7000, 6000,
    )
    anchor_angle = np.arctan(y_x[0]) * 180 / np.pi
    plt.plot(x, y, label=f"no slackliners, 6000N, anchor angle {anchor_angle:.2f}deg")

    x, y, n, y_x, n_x, T = integrate_natural_length(
        ideal, dyneemite_pro, 7000, n[-1],
        masses=[(10, 70)],
    )
    anchor_angle = np.arctan(y_x[0]) * 180 / np.pi
    plt.plot(x, y, label=f"1 slackliner, 70kg, 10m, {T}N, anchor angle {anchor_angle:.2f}deg")
    plt.legend()
    plt.show()


    # Now use the natural length n[-1] to add a slackliner of 70kg at 50m:
    #x, y, n, y_x, n_x, anchor_tension = integrate_natural_length(
    #    ideal, dyneemite_pro, 100, n[-1],
    #    masses=[(50, 70)],
    #)
    #plt.plot(x, y, label=f"1x70kg, {anchor_tension}N")
    
    # 9m gap, 9.7m natural, 70kg slacklinerat at 3m
