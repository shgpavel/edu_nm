import numpy as np
import pycuda.autoinit
import pycuda.driver as drv

from pycuda.compiler import SourceModule
from numba import njit

from aux.var import A, B, XI
from aux.func import func
from aux.init_h import init_h

S = 2

with open("methods/kernels/rk2.cu", "r") as f:
    kernel_code = f.read()

module = SourceModule(kernel_code, options=["-use_fast_math"])
rk2_cuda = module.get_function("rk2_cuda")

def rk2(y0, t0, t_end, h):
    if len(y0) % 2 != 0:
        raise ValueError("try something else")

    a21 = XI
    b2 = 1 / (2 * XI)
    b1 = 1 - b2
    
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
    
    rk2_cuda(
        np.float64(h),
        np.float64(A),
        np.float64(B),
        np.float64(a21),
        np.float64(b1),
        np.float64(b2),
        np.int32(num_steps),
        y_initial_gpu,
        y_values_gpu,
        block=(threads_per_block, 1, 1),
        grid=(blocks, 1)
    )
    
    drv.memcpy_dtoh(y_values_host, y_values_gpu)
    
    return t_values, y_values_host

@njit
def rk2_cpu(y0, t0, t_end, h):

    a21 = XI
    b2 = 1 / (2 * XI)
    b1 = 1 - b2

    num_steps = int(np.ceil((t_end - t0) / h)) + 1
    t_values = np.linspace(t0, t_end, num_steps)
    y_values = np.zeros((num_steps, len(y0)), dtype=np.float64)

    y_values[0, :] = y0

    k1 = np.zeros(len(y0), dtype=np.float64)
    k2 = np.zeros(len(y0), dtype=np.float64)
    y_temp = np.zeros(len(y0), dtype=np.float64)
    
    for i in range(1, num_steps):
        t = t_values[i - 1]
        y = y_values[i - 1]

        k1 = func(t, y)
        y_temp = y + a21 * h * k1
        k2 = func(t + a21 * h, y_temp)

        y_new = y + h * (b1 * k1 + b2 * k2)
        y_values[i] = y_new

    return t_values, y_values

@njit
def rk2_tol(y0, t0, t_end, tol):
    h = init_h(y0, t0, t_end, tol, S)
    sol, ri_hat = np.ones((2, len(y0)), dtype=np.float64)
    
    while np.any(np.abs(ri_hat) > tol):
        t, y = rk2_cpu(y0, t0, t_end, h)
        h *= 0.5
        t_half, y_half = rk2_cpu(y0, t0, t_end, h)

        ri_hat = (y_half[-1] - y[-1]) / (pow(2, S) - 1)
        sol = y_half[-1] + ri_hat

    return np.array([t_end]), sol[np.newaxis, :]
