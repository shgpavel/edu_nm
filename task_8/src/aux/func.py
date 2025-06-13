import numpy as np
import math

from aux.var import A, B


def func(t, y):
    return np.array([A * y[1], -B * y[0]], dtype=np.float64)


def analytical_solution(t):
    y1 = (7 * math.sqrt(35) * math.pi * math.sin(math.sqrt(35) * t / 12)) / 60 + (
        5 * math.pi * math.cos(math.sqrt(35) * t / 12)
    ) / 12

    y2 = (7 * math.pi * math.cos(math.sqrt(35) * t / 12)) / 12 - (
        5 * math.sqrt(35) * math.pi * math.sin(math.sqrt(35) * t / 12)
    ) / 84

    return np.array([y1, y2], dtype=np.float64)
