import tensorflow as tf
from findiff import FinDiff

def tensor(N, dx, order=1, accuracy=4, dtype=tf.float32):
    """
    Returns a tensor that calculates the derivative of a function using a
    high-precision finite difference stencil. Can be used to calculate
    derivatives in a way that tensorflow can automatically differentiate.

    Inputs:
        N: the number of points in discrete grid
        dx: the step size of the grid
        order: the order of the derivative, defaults to 1
        accuracy: the accuracy of the stencil, defaults to 4

    Returns:
        A tensor that calculates the derivative of a function.
    """
    dx_dt_scipy = FinDiff(0, dx, order, acc=accuracy).matrix((N,))
    return tf.constant(dx_dt_scipy.todense(), dtype=dtype)
