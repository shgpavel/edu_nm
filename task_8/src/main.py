import numpy as np
import time
import matplotlib.pyplot as plt

from aux.var import A, B
from methods.rk2 import rk2, rk2_cpu, rk2_tol
from methods.rk4 import rk4, rk4_cpu, rk4_tol
from methods.rk2wc import rk2wc_cpu
from methods.rk4wc import rk4wc_cpu

EPSILON = 1e-4
RHO = 1e-5
TOL_NAMES = ["rk2wc_cpu", "rk4wc_cpu", "rk2_tol", "rk4_tol"]

def print_perf(*funcs, y0, t0, t_end, h, ntests=20):
    for i, f in enumerate(funcs):
        execution_times = []
        for _ in range(ntests):
            start_time = time.time()
            
            if f.__name__ in TOL_NAMES:
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
        
        print(f"{f.__name__}\n"
              f"{'avg':<10}{execution_time / ntests:.6f} sec\n"
              f"{'low 1p':<10}{low_1p:.6f} sec\n"
              f"{'95p':<10}{high_95p:.6f} sec")

        if i != len(funcs) - 1:
            print("\n", end="")

def print_sols(*funcs, y0, t0, t_end, h=None): 
    for f in funcs:
        if f.__name__ in TOL_NAMES:
            t_values, y_values = f(y0, t0, t_end, RHO)    
        else:
            t_values, y_values = f(y0, t0, t_end, h)

        y1_pi, y2_pi = y_values[-1, 0], y_values[-1, 1]
        print(f"{f.__name__:<10} y1(pi) = {y1_pi:.6f}, y2(pi) = {y2_pi:.6f}")
        
def main():
    y0 = np.array([B * np.pi, A * np.pi], dtype=np.float64)
    t0 = 0.0
    t_end = np.pi
    h = 0.001
    
    print_sols(rk2, rk2_cpu, rk2wc_cpu, rk4, rk4_cpu, rk4wc_cpu,
               rk2_tol, rk4_tol, y0=y0, t0=t0, t_end=t_end, h=h)
    print("\n", end="")
    
    print_perf(rk2, rk2_cpu, rk2wc_cpu, rk4, rk4_cpu, rk4wc_cpu,
               rk2_tol, rk4_tol, y0=y0, t0=t0, t_end=t_end, h=h)
      
if __name__ == "__main__":
    main()
