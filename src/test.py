import matplotlib.pyplot as plt
import numpy as np
from .lagrangians import ideal, first_approximation
from .slacklines import dyneemite_pro
from .integrator import integrate

def main():
    # Run the integrator for an ideal dyneemite_pro slackline, for tensions
    # 1000, 1500, 2000, ..., 4000, and angle 0.1,0.2,0.3,0.4,0.5:
    for T0 in [1000, 2000, 3000]:
        for theta0 in [0.1, 0.2, 0.3]:
            print(f"Integrating for T0={T0}, theta0={theta0}")
            x, y, n, y_x, n_x = integrate(
                ideal,
                dyneemite_pro,
                T0,
                theta0,
                masses=[(5.0, 70.0), (10.0, 80)],
                length_cutoff=100000.0,
            )
            plt.plot(x, y, label=f"T0={T0}, theta0={theta0}")
    plt.legend()
    plt.show()
