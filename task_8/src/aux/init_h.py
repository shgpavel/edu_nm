import numpy as np

from aux.func import func
from numba import njit

H_TOL = 1e-2

@njit
def init_h(y0: np.float64, t0: np.float64, t_end: np.float64,
           tol: np.float64, s: np.int32) -> np.float64:

    delta = pow((1 / max(abs(t0), abs(t_end))), s + 1) + \
        pow(np.linalg.norm(func(t0, y0), ord=np.inf), s + 1)
    h = pow(tol / delta, 1/(s + 1))

    if np.linalg.norm(y0, ord=2) < H_TOL:
        y0 += h * func(t0, y0)
        t0 += h
        new_h = init_h(y0, t0, t_end, tol, s)
        return min(h, new_h)
    
    return h
