import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
from findiff import FinDiff

# Constants:
C = 10 # number of cells
G = 9.81 # m / s^2 (gravitational acceleration)
K = 2500.0 # N (newtons per 1% stretch in webbing)
M = 0.044 # kg / m (mass density of rope in natural units)
X = 100.0 # m (gap length)
N = 100.0 # m (natural length of webbing) 

# Create our C mass cells. Each mass cell has weight N*m/C and is connected to
# its neighbours with a segment of webbing with stretch factor K:
m = np.concatenate(([0], np.ones(C) * N * M / C, [0]))
x = np.linspace(0, X, C+2)
y = np.zeros(C+2)
x_t = np.zeros(C+2)
y_t = np.zeros(C+2)

# Create our C+1 webbing segments. The tension in each segment is given by
# T = K * (l - n) / n, where n is the natural length of the segment and l is
# the current length of the segment.
n = np.ones(C+1) * N / (C + 1)

# Book-keeping for efficiency:
l = np.zeros(C+1) # length of each webbing segment
t = np.zeros(C+1) # tension in each webbing segment
t_x = np.zeros(C+1) # x-component of tension in each webbing segment
t_y = np.zeros(C+1) # y-component of tension in each webbing segment
f_x = np.zeros(C+1) # x-component of force on each mass cell
f_y = np.zeros(C+1) # y-component of force on each mass cell


def iterate(dt):
    global x, y, x_t, y_t, l, t, t_x, t_y, f_x, f_y

    # Calculate the current length of each webbing segment:
    l = np.sqrt((x[1:] - x[:-1])**2 + (y[1:] - y[:-1])**2)

    # Calculate the tension in each webbing segment. We use a high-precision
    # finite difference method to stabilise the integration:
    t = 100.0 * K * (l - n) / n

    # Calculate the horizontal and vertical components of the tensions:
    t_x = t * (x[1:] - x[:-1]) / l
    t_y = t * (y[1:] - y[:-1]) / l

    # The total horizontal force on a mass cell is the difference between the
    # tensions on its left and right webbing segments:
    f_x = t_x[1:] - t_x[:-1]

    # The total vertical force on a mass cell is the difference between the
    # tensions on its left and right webbing segments, plus the weight of the
    # mass cell:
    f_y = t_y[1:] - t_y[:-1] - m[1:-1] * G

    # Calculate the velocity of each mass cell:
    x_t[1:-1] += dt * f_x / m[1:-1]
    y_t[1:-1] += dt * f_y / m[1:-1]

    # Calculate the position of each mass cell:
    x[1:-1] += dt * x_t[1:-1]
    y[1:-1] += dt * y_t[1:-1]

def main():
    # Create an animated plot of x, y over time, where the plot is updated
    # every 0.1 seconds:
    fig, ax = plt.subplots()
    ax.set_ylim(-100, 5)
    line, = ax.plot(x, y)
    def animate(i):
        iterate(0.01)
        line.set_data(x, y)
        return line,
    ani = animation.FuncAnimation(fig, animate, interval=100)
    plt.show()
