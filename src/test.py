import matplotlib.pyplot as plt
from .api import Constraints

def main():
    c = Constraints(gap_length=100, anchor_tension=1000)
    c.add_slackliner(position=50, mass=90)
    r = c.rig()
    r.plot()
    plt.show()
