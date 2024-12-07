import numpy as np
import pycuda.autoinit
import pycuda.driver as drv

from pycuda.compiler import SourceModule
from numba import njit

from aux.func import func
from aux.var import A, B


with open("methods/kernels/rk4.cu", "r") as f:
    kernel_code = f.read()

module = SourceModule(kernel_code, options=["-use_fast_math"])
rk4_cuda = module.get_function("rk4_cuda")

def rk4(y0, t0, t_end, h):
    if len(y0) % 2 != 0:
        raise ValueError("try something else")
    
    y0 = np.array(y0, dtype=np.float32)
    N = len(y0) // 2

    num_steps = int(np.ceil((t_end - t0) / h)) + 1
    t_values = np.linspace(t0, t_end, num_steps, dtype=np.float32)
    
    y_values_host = np.zeros((num_steps, 2 * N), dtype=np.float32)
    y_values_host[0, :] = y0

    y_initial_gpu = drv.mem_alloc(y0.nbytes)
    drv.memcpy_htod(y_initial_gpu, y0)

    y_values_gpu = drv.mem_alloc(y_values_host.nbytes)
    
    threads_per_block = 2
    blocks = (N + threads_per_block - 1) // threads_per_block
    
    rk4_cuda(
        np.float32(h),
        np.float32(A),
        np.float32(B),
        np.int32(num_steps),
        y_initial_gpu,
        y_values_gpu,
        block=(threads_per_block, 1, 1),
        grid=(blocks, 1)
    )
    
    drv.memcpy_dtoh(y_values_host, y_values_gpu)
    
    return t_values, y_values_host

@njit
def rk4_cpu(y0, t0, t_end, h):

    num_steps = int(np.ceil((t_end - t0) / h)) + 1
    t_values = np.linspace(t0, t_end, num_steps)
    y_values = np.zeros((num_steps, len(y0)), dtype=np.float32)

    y_values[0, :] = y0

    k1 = np.zeros(len(y0), dtype=np.float32)
    k2 = np.zeros(len(y0), dtype=np.float32)
    y_temp = np.zeros(len(y0), dtype=np.float32)
    
    for i in range(1, num_steps):
        t = t_values[i - 1]
        y = y_values[i - 1]

        k1 = h * func(t, y)
        k2 = h * func(t + 0.5 * h, y + 0.5 * k1)
        k3 = h * func(t + 0.5 * h, y + 0.5 * k2)
        k4 = h * func(t + h, y + k3)

        y_new = y + 1/6 * (k1 + 2 * k2 + 2 * k3 + k4)
        y_values[i] = y_new

    return t_values, y_values
