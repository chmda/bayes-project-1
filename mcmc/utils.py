import numpy as np
from scipy.linalg import cholesky


def multivariate_normal_p(mean, precision, size=None):
    """Draw random samples from a multivariate normal distribution with precision matrix.

    :param mean: Mean
    :param precision: Precision matrix
    :param size: Shape, defaults to None
    :return: The drawn samples
    """
    mean = np.asarray(mean)
    precision = np.asarray(precision)
    if size is None:
        shape = []
    elif isinstance(size, (int, np.integer)):
        shape = [size]
    else:
        shape = size

    if len(mean.shape) != 1:
        raise ValueError("mean must be 1 dimensional")
    if (len(precision.shape) != 2) or (precision.shape[0] != precision.shape[1]):
        raise ValueError("precision must be 2 dimensional and square")
    if mean.shape[0] != precision.shape[0]:
        raise ValueError("mean and precision must have same length")

    final_shape = list(shape[:])
    final_shape.append(mean.shape[0])
    x = np.random.standard_normal(final_shape).reshape(-1, mean.shape[0])
    U = cholesky(precision, lower=False)

    x = np.linalg.solve(U, x)
    x += mean
    x.shape = tuple(final_shape)
    return x
