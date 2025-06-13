import numpy as np

try:
    import pycuda.autoinit
    import pycuda.driver as drv

    from pycuda.compiler import SourceModule

    pycuda_avail = True
except ImportError:
    pycuda_avail = False

from aux.var import A, B
from aux.func import func
from aux.init_h import init_h

S = 4

if pycuda_avail:
    with open("methods/kernels/rk4.cu", "r") as f:
        kernel_code = f.read()

    module = SourceModule(kernel_code, options=["-use_fast_math"])
    rk4_cuda = module.get_function("rk4_cuda")


def rk4(y0, t0, t_end, h):
    if len(y0) % 2 != 0:
        raise ValueError("This function only doing systems of 2")

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
        grid=(blocks, 1),
    )

    drv.memcpy_dtoh(y_values_host, y_values_gpu)

    return t_values, y_values_host


def rk4_cpu(y0, t0, t_end, h):
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
        k2 = func(t + 0.5 * h, y + 0.5 * h * k1)
        k3 = func(t + 0.5 * h, y + 0.5 * h * k2)
        k4 = func(t + h, y + h * k3)

        y_new = y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        y_values[i] = y_new

    return t_values, y_values


def rk4_tol(y0, t0, t_end, tol, flag=0):
    h = init_h(y0, t0, t_end, tol, S)
    sol, ri_hat = np.ones((2, len(y0)), dtype=np.float64)

    i = 0  # i fixes ri_hat
    while np.any(np.abs(ri_hat) > tol) or i < 5:
        t, y = rk4_cpu(y0, t0, t_end, h)
        h *= 0.5
        t_half, y_half = rk4_cpu(y0, t0, t_end, h)

        ri_hat = (y_half[-1] - y[-1]) / (pow(2, S) - 1)
        sol = y_half[-1] + ri_hat
        i += 1

    if flag == 0:
        return np.array([t_end]), sol[np.newaxis, :]
    return h, t_half, y_half
