import matplotlib.pyplot as plt
import tensorflow as tf
import numpy as np
from findiff import FinDiff

# Copilot: use float64!

# Hyperparameters:
M = 1000 # number of training steps
N = 10 # number of time steps
optimiser = tf.keras.optimizers.Adam(learning_rate=0.1)

# We're going to run a test. First, we'll create a random set of 'x' values,
# and use the optimiser to try and minimise it's time derivative. If this works,
# we should end up with a set of 'x' values that are all the same.
dt = 0.1
x = tf.Variable(tf.random.uniform((N, 1), dtype=tf.float64))
t = tf.cumsum(tf.ones((N, 1), dtype=tf.float64)*dt, axis=0)

# We use findiff for computing derivatives, but to play nicely with tensorflow
# we have to convert the stencil into a tensor:
dx_dt_scipy = FinDiff(0, dt, acc=2).matrix((N,))
dx_dt = tf.constant(dx_dt_scipy.todense())

# We'll use a simple loss function that just sums the squares of the time
# derivative of 'x'.
def loss():
    dx = tf.matmul(dx_dt, x)
    return tf.reduce_sum(dx**2)

# Now we run the training loop:
for i in range(M):
    optimiser.minimize(loss, [x])
    if i % 100 == 0:
        print(f"Step {i}: {x.numpy()}")

# Check that we got the right answer:
print(f"Final answer: {x.numpy()}")
