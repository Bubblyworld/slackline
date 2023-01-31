import matplotlib.pyplot as plt
import numpy as np
from .calculus import first_order_euler_lagrange, mass_boundary_conditions
from .lagrangians import ideal, first_approximation
from .slacklines import dyneemite_pro
from .integrator import integrate, integrate_length_tension, integrate_natural_length

def main():
    # Precomputation for speed:
    print("Precomputing...")
    lagrangian = ideal
    slackline = dyneemite_pro
    foel = first_order_euler_lagrange(lagrangian)
    masses = [(50.0, 70.0), (100.0, 70.0)] # (x, M)
    mbcs = [mass_boundary_conditions(lagrangian, x, slackline.m,
                                     slackline.g, slackline.K, M)
                for x, M in masses]

    # Solve for the length-tension curve:
    print("Solving for length-tension curves...")
    for T0 in [2000, 2500, 3000, 3500, 4000]:
        print(f"Integrating for T0={T0}")
        x, y, n, y_x, n_x = integrate_length_tension(
            lagrangian,
            slackline,
            130,
            T0,
            masses=masses,
            foel=foel,
            mbcs=mbcs,
        )
        plt.plot(x, y, label=f"T0={T0}")
    plt.legend()
    plt.show()
