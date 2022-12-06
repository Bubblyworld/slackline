import matplotlib.pyplot as plt
import numpy as np
from .lagrangians import ideal, first_approximation
from .slacklines import dyneemite_pro
from .integrator import integrate, integrate_length_tension

def main():
    x, y, n, y_x, n_x = integrate_length_tension(ideal, dyneemite_pro, 100.0, 1000.0)
    plt.plot(x, y)
    plt.show()
