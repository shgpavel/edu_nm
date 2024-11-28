import numpy as np
import pycuda.autoinit
import pycuda.driver as drv

from aux.var import A, B, XI
from pycuda.compiler import SourceModule

S = 2

mod = SourceModule("""
#include <math.h>
__global__ void rk2_step(
        float h,
        unsigned num_steps,
        float A,
        float B,
        float a21,
        float b1,
        float b2,
        float *y_initial,
        float *y_values) {
    
    int system_idx = blockIdx.x * blockDim.x + threadIdx.x;

    float *y_current = y_initial + 2 * system_idx;
    float *y_all_steps = y_values + system_idx;

    float y1 = y_current[0];
    float y2 = y_current[1];

    y_all_steps[0] = y1;
    y_all_steps[num_steps * 2 * gridDim.x] = y2;

    for (unsigned step = 1; step < num_steps; ++step) {
        float k1_y1 = A * y2;
        float k1_y2 = -B * y1;
            
        float y_temp1 = y1 + a21 * h * k1_y1;
        float y_temp2 = y2 + a21 * h * k1_y2;
            
        float k2_y1 = A * y_temp2;
        float k2_y2 = -B * y_temp1;
            
        y1 = y1 + h * (b1 * k1_y1 + b2 * k2_y1);
        y2 = y2 + h * (b1 * k1_y2 + b2 * k2_y2);
            
        y_all_steps[step * 2 * gridDim.x] = y1;
        y_all_steps[step * 2 * gridDim.x + 1] = y2;
    }
}""", options=["-use_fast_math"])


rk2_step = mod.get_function("rk2_step")

def rk2(y0, t0, t_end, h):
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
    
    threads_per_block = 64
    blocks = (N + threads_per_block - 1) // threads_per_block
    
    rk2_step(
        np.float32(h),
        np.int32(num_steps),
        np.float32(A),
        np.float32(B),
        np.float32(XI),
        np.float32(1.0 - 1.0/(2.0 * XI)),
        np.float32(1.0/(2.0 * XI)),
        y_initial_gpu,
        y_values_gpu,
        block=(threads_per_block, 1, 1),
        grid=(blocks, 1)
    )
    
    drv.memcpy_dtoh(y_values_host, y_values_gpu)
    
    return t_values, y_values_host


def rk2wc(f, y0, t0, t_end, tol):
    delta = pow((1 / max(abs(t0), abs(t_end))), S)
    # + pow(np.linalg.norm(f, 2), S + 1)
    h = pow(tol / delta, 1/(S + 1))
    
    a21 = XI
    b2 = 1 / (2 * XI)
    b1 = 1 - b2

    t_values = np.arange(t0, t_end + h, h)
    y_lin_values = np.zeros((len(t_values), len(y0)))
    y_cap_values = np.zeros((len(t_values), len(y0)))
    # y_values = np.zeros((len(t_values), len(y0)))
    
    y_lin_values[0, :], y_cap_values[0, :] = y0, y0
    
    for i in range(1, len(t_values)):
        t = t_values[i - 1]
        y = y_values[i - 1]
        
        k1 = h * f(t, y)
        k2 = h * f(t + XI * h, y + a21 * k1 * h)
        
        y_values[i, :] = y + b1 * k1 + b2 * k2

    return t_values, y_values
