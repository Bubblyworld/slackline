import matplotlib.pyplot as plt
import numpy as np
from .calculus import first_order_euler_lagrange, mass_boundary_conditions
from .lagrangians import ideal, first_approximation
from .slacklines import dyneemite_pro
from .integrator import integrate, integrate_length_tension

def main():
    x, y, n, y_x, n_x = integrate_length_tension(
        ideal, dyneemite_pro, 100.0, 3500.0,
        masses=[(70.0, 80.0)],
    )
    plt.plot(x, y)
    plt.show()
