import numpy as np
import pycuda.autoinit
import pycuda.driver as drv

from aux.var import A, B, XI
from pycuda.compiler import SourceModule

S = 2

mod = SourceModule("""
#include <math.h>

__global__ void rk2wc_cuda(
        float h,
        float A,
        float B,
        float a21,
        float b1,
        float b2,
        float tol,
        unsigned num_steps,
        float *y_initial,
        float *y_values) {
    
    int system_idx = blockIdx.x * blockDim.x + threadIdx.x;

    float *y_current = y_initial + 2 * system_idx;
    float *y_all_steps = y_values + system_idx;

    float y1 = y_current[0];
    float y2 = y_current[1];

    y_all_steps[0] = y1;
    y_all_steps[num_steps * 2 * gridDim.x] = y2;

    float h_2 = h * 0.5f;

    for (unsigned step = 1; step < num_steps; ++step) {
        /* f */
        float f_y1 = A * y2;
        float f_y2 = -B * y1;

        /* associative part */
        float ap_1 = a21 * h * f_y1;
        float ap_2 = a21 * h * f_y2;

        /* first part of full step */
        float y_1p1 = y1 + ap_1;
        float y_2p1 = y2 + ap_2;

        /* apply f 2nd time */
        float f2_y1 = A * y_1p1;
        float f2_y2 = -B * y_2p1;
        
        /* first part of full step */
        float y_1p1_2 = y1 + ap_1 / 2;
        float y_2p1_2 = y2 + ap_2 / 2;

        /* apply f 2nd time */
        float f2_y1_2 = A * y_1p1_2;
        float f2_y2_2 = -B * y_2p1_2;

        float y1_full = y1 + b1 * h * f_y1 + b2 * h * f2_y1;
        float y2_full = y2 + b1 * h * f_y2 + b2 * h * f2_y2;
        
        float y1_half = y1 + b1 * h_2 * f_y1 + b2 * h_2 * f2_y1_2;
        float y2_half = y2 + b1 * h_2 * f_y1 + b2 * h_2 * f2_y2_2;

        float error1 = fabsf(y1_full - y1_half);
        float error2 = fabsf(y2_full - y2_half);
        float error = fmaxf(error1, error2) / 3.0f;

        y_all_steps[step * 2 * gridDim.x] = y1;
        y_all_steps[step * 2 * gridDim.x + 1] = y2; 

        if (error > tol) {
            h *= 0.5f;
        } else {
            
        }
    }
}
""", options=["-use_fast_math"])

 rk2wc_cuda = mod.get_function("rk2wc_cuda")    

def rk2wc(f, y0, t0, t_end, tol):
    if len(y0) % 2 != 0:
        raise ValueError("try something else")
    
    delta = pow((1 / max(abs(t0), abs(t_end))), S)
    # + pow(np.linalg.norm(f, 2), S + 1)
    h = pow(tol / delta, 1/(S + 1))


    y0 = np.array(y0, dtype=np.float32)
    N = len(y0) // 2
    
    a21 = XI
    b2 = 1 / (2 * XI)
    b1 = 1 - b2
    
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
        np.float32(XI),
        np.float32(1.0 - 1.0/(2.0 * XI)),
        np.float32(1.0/(2.0 * XI)),
        np.int32(num_steps),
        y_initial_gpu,
        y_values_gpu,
        block=(threads_per_block, 1, 1),
        grid=(blocks, 1)
    )
    
    drv.memcpy_dtoh(y_values_host, y_values_gpu)
    return t_values, y_values_host
