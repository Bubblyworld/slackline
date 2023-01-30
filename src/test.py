import matplotlib.pyplot as plt
import numpy as np
from .calculus import first_order_euler_lagrange, mass_boundary_conditions
from .lagrangians import ideal, first_approximation
from .slacklines import dyneemite_pro
from .integrator import integrate, integrate_length_tension, integrate_natural_length

def main():
    x, y, n, a, b = integrate_length_tension(ideal, dyneemite_pro, 100, 4000, masses=[(50, 80)])
    plt.plot(x, y)
    plt.show()
