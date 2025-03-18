import numpy as np
import time
import matplotlib.pyplot as plt

try:
    import pycuda.autoinit
    pycuda_avail = True
except ImportError:
    pycuda_avail = False

from aux.var import A, B
from aux.func import analytical_solution
from methods.rk2 import rk2, rk2_cpu, rk2_tol
from methods.rk4 import rk4, rk4_cpu, rk4_tol
from methods.rk2wc import rk2wc_cpu
from methods.rk4wc import rk4wc_cpu

EPSILON = 1e-4
RHO = 1e-5
TOL_NAMES = ["rk2_tol", "rk4_tol"]
COR_NAMES = ["rk2wc_cpu", "rk4wc_cpu"]
CUDADEP_NAMES = ["rk2", "rk4"]

def print_perf(*funcs, y0, t0, t_end, h, ntests=100):
    for i, f in enumerate(funcs):
        execution_times = []

        if f.__name__ in CUDADEP_NAMES and not(pycuda_avail):
            continue
        elif f.__name__ in COR_NAMES:
            last_arg = RHO
        elif f.__name__ in TOL_NAMES:
            last_arg = EPSILON
        else:
            last_arg = h
        
        for _ in range(ntests):
            start_time = time.time()

            f(y0, t0, t_end, last_arg)
            
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

def print_sols(*funcs, y0, t0, t_end, h=None, actsol): 
    for f in funcs:
        if f.__name__ in CUDADEP_NAMES and not(pycuda_avail):
            continue
        elif f.__name__ in COR_NAMES:
            t_values, y_values = f(y0, t0, t_end, RHO)
        elif f.__name__ in TOL_NAMES:
            t_values, y_values = f(y0, t0, t_end, EPSILON)
        else:
            t_values, y_values = f(y0, t0, t_end, h)

        y1_pi, y2_pi = y_values[-1, 0], y_values[-1, 1]
        diff = np.linalg.norm(actsol - [y1_pi, y2_pi], ord=2)
        print(f"{f.__name__:<10} y1(pi) = {y1_pi:.6f}", end="")
        print(f" y2(pi) = {y2_pi:.6f} diff = {diff:e}")

def graph_ftol(y0, t0, t_end):
    hopt2, t2, y2 = rk2_tol(y0, t0, t_end, EPSILON, 1)
    hopt4, t4, y4 = rk4_tol(y0, t0, t_end, EPSILON, 1)
    print(f"optimal step rk2 {hopt2:e}\noptimal step rk4 {hopt4:e}")

    plt.figure(figsize=(10, 6))

    t = np.linspace(0, np.pi, 10000)
    y = np.array([analytical_solution(ti) for ti in t])

    plt.plot(t, y, 'b-.', label=f"y")
    plt.plot(t2, y2, label=f'y_rk2')
    plt.plot(t4, y4, label=f'y_rk4')
    
    plt.xlabel('t')
    plt.ylabel('y')
    plt.title('global')
    plt.legend()
    plt.grid(True)
    plt.show()

def graph_ht(y0, t0, t_end):
    t2, h2 = rk2wc_cpu(y0, t0, t_end, RHO, 1)
    t4, h4 = rk4wc_cpu(y0, t0, t_end, RHO, 1)
    
    plt.figure(figsize=(10, 6))
    plt.plot(t2, h2, label=f'h2')
    plt.plot(t4, h4, label=f'h4')

    plt.xlabel('t')
    plt.ylabel('h')
    plt.title('ht')
    plt.legend()
    plt.grid(True)
    plt.show()

def graph_cnt(y0, t0, t_end):
    tol = 1e-1
    rk2wc_cnt, rk4wc_cnt = [], []
    while tol >= 1e-13:
        rk2wc_cntd = rk2wc_cpu(y0, t0, t_end, tol, 2)
        rk4wc_cntd = rk4wc_cpu(y0, t0, t_end, tol, 2)
        
        rk2wc_cnt.append(rk2wc_cntd)
        rk4wc_cnt.append(rk4wc_cntd)
        
        tol /= 10

    x = [i + 1 for i in range(len(rk2wc_cnt))]
    plt.figure(figsize=(10, 6))

    plt.plot(x, rk2wc_cnt, label=f"rk2wc")
    plt.plot(x, rk4wc_cnt, label=f"rk4wc")
    
    plt.xlabel('drp')
    plt.ylabel('cnt')
    plt.title('cnt')
    plt.legend()
    plt.grid(True)
    plt.show()

def graph_loc(y0, t0, t_end):
    plt.figure(figsize=(10, 6))
    t_end *= 0.3
    
    t = np.linspace(0, t_end, 10000)
    y = np.array([analytical_solution(ti) for ti in t])
    t1, y1 = rk2wc_cpu(y0, t0, t_end, RHO)
    ta, ha, ea = rk2wc_cpu(y0, t0, t_end, RHO, 4)

    plt.plot(t, y, label=f"y")
    plt.plot(t1, y1, label=f"y1")

    for i in range(6):
        ti, yi = rk2wc_cpu(y0, t0 + ha[i], t_end, RHO)
        plt.plot(ti, yi, label=f"y{i}")
        
    plt.xlabel('t')
    plt.ylabel('y')
    plt.title('loc')
    plt.legend()
    plt.grid(True)
    plt.show()
    
def main():
    y0 = np.array([B * np.pi, A * np.pi], dtype=np.float64)
    t0 = 0.0
    t_end = np.pi
    h = 1e-4
    
    actsol = analytical_solution(np.pi) 
    print(f"{"an sol":<10} y1(pi) = {actsol[0]:.6f} y2(pi) = {actsol[1]:.6f}")
    print_sols(rk2, rk2_cpu, rk2wc_cpu, rk4, rk4_cpu, rk4wc_cpu,
               rk2_tol, rk4_tol, y0=y0, t0=t0, t_end=t_end, h=h, actsol=actsol)
    print("\n", end="")
 

    # rk2_cpu, rk4_cpu,
    print_perf(rk2, rk2wc_cpu, rk4, rk4wc_cpu,
               rk2_tol, rk4_tol, y0=y0, t0=t0, t_end=t_end, h=h)
    print("\n", end="")
    
    '''
    graph_ftol(y0, t0, t_end)
    graph_ht(y0, t0, t_end)
    
    graph_loc(y0, t0, t_end)
    graph_cnt(y0, t0, t_end)
    '''
if __name__ == "__main__":
    main()
