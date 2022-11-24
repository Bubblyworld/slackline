import tensorflow as tf
from src.derivative import tensor

# We run some simple tests on a discrete grid of 100 points:
N = 10
dt = 0.1
d_dt = tensor(N, dt, order=1, accuracy=4, dtype=tf.float64)
dd_dt = tensor(N, dt, order=2, accuracy=4, dtype=tf.float64)
threshold = 1e-4

def test_constant():
    """
    Tests that the derivative of a constant function is zero.
    """
    x = tf.ones((N, 1), dtype=tf.float64)
    dx_dt = tf.matmul(d_dt, x)
    assert tf.reduce_sum(tf.abs(dx_dt))/N < threshold

def test_linear():
    """
    Tests that the derivative of a linear function is a constant, and that
    the second derivative is zero.
    """
    x = tf.linspace(0.0, 1.0, N)[:, None]
    x = tf.cast(x, dtype=tf.float64)
    dx_dt = tf.matmul(d_dt, x)
    assert tf.reduce_sum(tf.abs(dx_dt - dx_dt[0]))/N < threshold
    ddx_dt = tf.matmul(dd_dt, x)
    assert tf.reduce_sum(tf.abs(ddx_dt))/N < threshold

def test_optimiser_constant():
    """
    Tests that minimising the first derivative of a random set of values gives
    a constant function. This makes sure tensorflow is playing nicely with the
    derivative tensor.
    """
    x = tf.Variable(tf.random.uniform((N, 1), dtype=tf.float64))
    optimiser = tf.keras.optimizers.Adam(learning_rate=0.1)
    for i in range(1000):
        optimiser.minimize(lambda: tf.reduce_sum(tf.matmul(d_dt, x)**2), [x])
    assert tf.reduce_sum(tf.abs(x - x[0]))/N < threshold
