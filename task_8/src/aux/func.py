import numpy as np

from numba import njit

from aux.var import A, B

@njit
def func(t, y):
    return np.array([A * y[1], -B * y[0]])
