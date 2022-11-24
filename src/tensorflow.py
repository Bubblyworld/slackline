import matplotlib.pyplot as plt
import tensorflow as tf
import numpy as np
from findiff import FinDiff

# Hyperparameters:
M = 1000 # number of training steps
N = 1000 # number of time steps
optimiser = tf.keras.optimizers.Adam(learning_rate=0.1)

# Constants:
m = 1.0 # kg (mass of the object)
g = 9.81 # m / s^2 (gravitational acceleration)
dt = 0.001 # s (time differential)
v0 = 0.0 # m (initial velocity)
x0 = 10.0 # m (initial position)
infinity = 1e6
