import numpy as np
import timeit
import matplotlib.pyplot as plt

from aux.var import A, B
from methods.rk2 import rk2
from methods.rk2wc import rk2wc

EPSILON = 1e-4

def main():
    y0 = [B * np.pi, A * np.pi]
    t0 = 0.0
    t_end = np.pi
    h = 0.0001

    t_values, y_values = rk2(y0, t0, t_end, h)
    y1_pi, y2_pi = y_values[-1, 0], y_values[-1, 1]
    print(f"y1(pi) = {y1_pi:.4f}, y2(pi) = {y2_pi:.4f}")

    execution_time = timeit.timeit(lambda: rk2(y0, t0, t_end, h),
                                   number=5)
    print(f"AVG: {execution_time / 100:.6f} sec")

    def f(t, y):
        return np.array([A * y[1], -B * y[0]])
    
    t_values, y_values = rk2wc(f, y0, t0, t_end, EPSILON)
    y1_pi, y2_pi = y_values[-1, 0], y_values[-1, 1]
    print(f"y1(pi) = {y1_pi:.4f}, y2(pi) = {y2_pi:.4f}")
    execution_time = timeit.timeit(lambda:
                                   rk2wc(f, y0, t0, t_end, EPSILON),
                                   number=5)
    print(f"AVG: {execution_time / 100:.6f} sec")


if __name__ == "__main__":
    main()
