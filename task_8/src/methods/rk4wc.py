import numpy as np

from aux.var import A, B, XI
from aux.func import func
from aux.init_h import init_h

S = 4
def rk4wc_cpu(y0, t0, t_end, tol, flag=0):
    t = t0
    y = y0

    t_values = [t]
    y_values = [y]
    
    h = init_h(y0, t0, t_end, tol, S)

    if flag == 1:
        ta, ha = [], []

    cnt = 0    
    while t < t_end:
        if t + h > t_end:
            h = t_end - t

        k1 = func(t, y)
        k2 = func(t + h/2, y + h/2 * k1)
        k3 = func(t + h/2, y + h/2 * k2)
        k4 = func(t + h, y + h * k3)
        y1 = y + (h/6) * (k1 + 2*k2 + 2*k3 + k4)

        h_half = h / 2

        k1 = func(t, y)
        k2 = func(t + h_half/2, y + h_half/2 * k1)
        k3 = func(t + h_half/2, y + h_half/2 * k2)
        k4 = func(t + h_half, y + h_half * k3)
        y_half = y + (h_half/6) * (k1 + 2*k2 + 2*k3 + k4)

        t += h_half        
        k1 = func(t, y_half)
        k2 = func(t + h_half/2, y_half + h_half/2 * k1)
        k3 = func(t + h_half/2, y_half + h_half/2 * k2)
        k4 = func(t + h_half, y_half + h_half * k3)
        y2 = y_half + (h_half/6) * (k1 + 2*k2 + 2*k3 + k4)

        cnt += 11
        error = np.linalg.norm((y2 - y1) / (1 - pow(2, -S)), ord=np.inf)
        eps2s = tol * pow(2, S)
        eps2s1 = tol / pow(2, S + 1)

        if error > eps2s:
            h = h_half
            t -= h_half
            
        elif error > tol:
            h = h_half
            y = y2
            t_values.append(t)
            y_values.append(y)
            
        elif error >= eps2s1:
            t += h_half
            y = y1
            t_values.append(t)
            y_values.append(y)

        else:
            h *= 2
            t += h_half
            y = y1
            t_values.append(t)
            y_values.append(y)

        if flag == 1:
            ta.append(t)
            ha.append(h)
            
    t_values = np.array(t_values)
    y_values = np.array(y_values, dtype=np.float64)

    if flag == 1:
        return ta, ha
    if flag == 2:
        return cnt
    return t_values, y_values
