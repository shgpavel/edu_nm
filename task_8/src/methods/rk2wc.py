import numpy as np
import pycuda.autoinit
import pycuda.driver as drv

from pycuda.compiler import SourceModule
from numba import njit

from aux.var import A, B, XI
from aux.func import func
from aux.init_h import init_h

S = 2

'''
with open("methods/kernels/rk2wc.cu", "r") as f:
    kernel_code = f.read()

module = SourceModule(kernel_code, options=["-use_fast_math"])
rk2wc_cuda = module.get_function("rk2wc_cuda")    

def rk2wc(y0, t0, t_end, tol):
    if len(y0) % 2 != 0:
        raise ValueError("try something else")

    a21 = XI
    b2 = 1 / (2 * XI)
    b1 = 1 - b2
    
    h = init_h(y0, t0, t_end, tol)
    
    y0 = np.array(y0, dtype=np.float64)
    N = len(y0) // 2
    
    num_steps = int(np.ceil((t_end - t0) / h)) + 1
    t_values = np.linspace(t0, t_end, num_steps, dtype=np.float64)
    
    y_values_host = np.zeros((num_steps, 2 * N), dtype=np.float64)
    y_values_host[0, :] = y0

    y_initial_gpu = drv.mem_alloc(y0.nbytes)
    drv.memcpy_htod(y_initial_gpu, y0)

    y_values_gpu = drv.mem_alloc(y_values_host.nbytes)
    
    threads_per_block = 64
    blocks = (N + threads_per_block - 1) // threads_per_block

    rk2wc_cuda(
        np.float64(t0),
        np.float64(t_end),
        np.float64(h),
        np.float64(tol),
        np.float64(A),
        np.float64(B),
        np.float64(a21),
        np.float64(b1),
        np.float64(b2),
        y_initial_gpu,
        y_values_gpu,
        block=(threads_per_block, 1, 1),
        grid=(blocks, 1)
    )
    
    drv.memcpy_dtoh(y_values_host, y_values_gpu)
    return t_values, y_values_host
'''

def rk2wc_cpu(y0, t0, t_end, tol):

    a21 = XI
    b2 = 1 / (2 * XI)
    b1 = 1 - b2
    
    t = t0
    y = np.copy(y0)
    t_values = [t]
    y_values = [np.copy(y)]

    h = init_h(y0, t0, t_end, tol, S)

    while t <= t_end:
        if t + h >= t_end:
            h = t_end - t

        k1 = func(t, y)
        y_temp = y + a21 * h * k1
        k2 = func(t + a21 * h, y_temp)
        y1 = y + h * (b1 * k1 + b2 * k2)
        
        h_half = h * 0.5

        y_temp = y + a21 * h_half * k1
        k2 = func(t + a21 * h_half, y_temp)
        y_half = y + h_half * (b1 * k1 + b2 * k2)

        y_temp = y_half + a21 * h_half * k1
        k2 = func(t + h_half + a21 * h_half, y_temp)
        y2 = y_half + h_half * (b1 * k1 + b2 * k2)
        
        error = np.linalg.norm(y2 - y1, ord=np.inf)    
        eps2s = tol * pow(2, S)
        eps2s1 = tol / pow(2, S + 1)

        if error > eps2s:
            h *= 0.5

        elif error > tol:
            t += h_half
            y = y2
            t_values.append(t)
            y_values.append(np.copy(y))
           
        elif error >= eps2s1:
            t += h
            y = y1
            t_values.append(t)
            y_values.append(np.copy(y))
            
        else:
            t += 2 * h
            y = y1
            t_values.append(t)
            y_values.append(np.copy(y))
            

    t_values = np.array(t_values)
    y_values = np.array(y_values, dtype=np.float64)
    
    return t_values, y_values
