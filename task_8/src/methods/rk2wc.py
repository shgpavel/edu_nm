import numpy as np
import pycuda.autoinit
import pycuda.driver as drv

from aux.var import A, B, XI
from aux.func import func

from pycuda.compiler import SourceModule


S = 2

with open("methods/kernels/rk2wc.cu", "r") as f:
    kernel_code = f.read()

module = SourceModule(kernel_code, options=["-use_fast_math"])
rk2wc_cuda = module.get_function("rk2wc_cuda")    

def rk2wc(f, y0, t0, t_end, tol):
    if len(y0) % 2 != 0:
        raise ValueError("try something else")

    a21 = XI
    b2 = 1 / (2 * XI)
    b1 = 1 - b2
    
    delta = pow((1 / max(abs(t0), abs(t_end))), S)
    # + pow(np.linalg.norm(f, 2), S + 1)
    h = pow(tol / delta, 1/(S + 1))

    y0 = np.array(y0, dtype=np.float32)
    N = len(y0) // 2
    
    num_steps = int(np.ceil((t_end - t0) / h)) + 1
    t_values = np.linspace(t0, t_end, num_steps, dtype=np.float32)
    
    y_values_host = np.zeros((num_steps, 2 * N), dtype=np.float32)
    y_values_host[0, :] = y0

    y_initial_gpu = drv.mem_alloc(y0.nbytes)
    drv.memcpy_htod(y_initial_gpu, y0)

    y_values_gpu = drv.mem_alloc(y_values_host.nbytes)
    
    threads_per_block = 64
    blocks = (N + threads_per_block - 1) // threads_per_block

    rk2wc_cuda(
        np.float32(h),
        np.float32(A),
        np.float32(B),
        np.float32(a21),
        np.float32(b1),
        np.float32(b2),
        np.int32(num_steps),
        y_initial_gpu,
        y_values_gpu,
        block=(threads_per_block, 1, 1),
        grid=(blocks, 1)
    )
    
    drv.memcpy_dtoh(y_values_host, y_values_gpu)
    return t_values, y_values_host


def rk2wc_cpu(y0, t0, t_end, tol=1e-5):
    
    a21 = XI
    b2 = 1 / (2 * XI)
    b1 = 1 - b2

    max_step_factor = 2
    min_step_factor = 0.1
    h_max = 1.0
    h_min = 1e-7
    
    t = t0
    y = np.array(y0, dtype=np.float32)
    t_values = [t]
    y_values = [y.copy()]

    def _rk_step(t, y, h):
        k1 = func(t, y)
        y_temp = y + a21 * h * k1
        k2 = func(t + a21 * h, y_temp)
        y_new = y + h * (b1 * k1 + b2 * k2)
        return y_new

    #def _init_h(t0, t_end, y0, tol):
    
    delta = pow((1 / max(abs(t0), abs(t_end))), S + 1) + \
        pow(np.linalg.norm(func(t, y), ord=np.inf), S + 1)
    h = pow(tol / delta, 1/(S + 1))
    
    while t < t_end:
        if t + h > t_end:
            h = t_end - t

        y1 = _rk_step(t, y, h)

        h_half = h * 0.5
        y_half = _rk_step(t, y, h_half)
        y2 = _rk_step(t + h_half, y_half, h_half)

        error_estimate = np.linalg.norm(y2 - y1, ord=np.inf)    
        scale = tol * max(np.linalg.norm(y2, ord=np.inf), 1.0)
        
        if error_estimate <= scale:
            t += h
            y = y2
            t_values.append(t)
            y_values.append(y.copy())
            
            if error_estimate < 1e-10:
                factor = max_step_factor
            else:
                factor = (scale / error_estimate) ** (1 / (S + 1))
                factor = min(max(factor, min_step_factor), max_step_factor)

            h = min(h * factor, h_max)
            continue

        factor = (scale / error_estimate) ** (1 / (S + 1))
        factor = max(factor, min_step_factor)
        h = max(h * factor, h_min)

        if h == h_min:
            t += h
            y = y2
            t_values.append(t)
            y_values.append(y.copy())

    t_values = np.array(t_values)
    y_values = np.array(y_values, dtype=np.float32)
    
    return t_values, y_values
