import numpy as np
import time
import matplotlib.pyplot as plt

from aux.var import A, B
from methods.rk2 import rk2, rk2_cpu, rk2_tol
from methods.rk2wc import rk2wc_cpu

EPSILON = 1e-4
RHO = 1e-5

def print_perf(f, y0, t0, t_end, h, ntests=100):
    execution_times = []
    for _ in range(ntests):
        start_time = time.time()

        if f.__name__ == "rk2wc_cpu":
            f(y0, t0, t_end, RHO)
        else:
            f(y0, t0, t_end, h)

        end_time = time.time()
        execution_times.append(end_time - start_time)
        
    execution_time = sum(execution_times)
    
    execution_times.sort()
    low_1p_count = max(1, ntests // 100)
    low_1p = np.mean(execution_times[-low_1p_count:])
    
    high_95p = np.percentile(execution_times, 95)
    low_1p2 = np.percentile(execution_times, 1)
    
    print(f"{f.__name__}\n"
          f"{'avg':<10}{execution_time / ntests:.6f} sec\n"
          f"{'low 1p':<10}{low_1p:.6f} sec\n"
          f"{'95p':<10}{high_95p:.6f} sec\n")

def main():
    y0 = np.array([B * np.pi, A * np.pi], dtype=np.float32)
    t0 = 0.0
    t_end = np.pi
    h = 0.0001

    t_values, y_values = rk2(y0, t0, t_end, h)
    y1_pi, y2_pi = y_values[-1, 0], y_values[-1, 1]
    print(f"{'rk2':<7} y1(pi) = {y1_pi:.6f}, y2(pi) = {y2_pi:.6f}")

    t_values, y_values = rk2_cpu(y0, t0, t_end, h)
    y1_pi, y2_pi = y_values[-1, 0], y_values[-1, 1]
    print(f"{'rk2_cpu':<7} y1(pi) = {y1_pi:.6f}, y2(pi) = {y2_pi:.6f}")

    t_values, y_values = rk2wc_cpu(y0, t0, t_end, RHO)
    y1_pi, y2_pi = y_values[-1, 0], y_values[-1, 1]
    print(f"{'rk2wc_cpu':<7} y1(pi) = {y1_pi:.6f}, y2(pi) = {y2_pi:.6f}\n")

    sol = rk2_tol(y0, t0, t_end, EPSILON)
    print(sol)
    #y1_pi, y2_pi = y_values[-1, 0], y_values[-1, 1]
    #print(f"{'rk2wc_cpu':<7} y1(pi) = {y1_pi:.6f}, y2(pi) = {y2_pi:.6f}\n")

    
    print_perf(rk2, y0, t0, t_end, h)
    print_perf(rk2_cpu, y0, t0, t_end, h)
    print_perf(rk2wc_cpu, y0, t0, t_end, h)
    
if __name__ == "__main__":
    main()
