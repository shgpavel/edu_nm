import numpy as np

from aux.var import A, B, XI
from aux.func import func
from aux.init_h import init_h

S = 4

def rk4wc_cpu(y0, t0, t_end, tol):

    t = t0
    y = np.copy(y0)
    t_values = [t]
    y_values = [np.copy(y)]

    h = init_h(y0, t0, t_end, tol, S)
    
    while t <= t_end:
        if t + h >= t_end:
            h = t_end - t

        k1 = h * func(t, y)
        k2 = h * func(t + 0.5 * h, y + 0.5 * k1)
        k3 = h * func(t + 0.5 * h, y + 0.5 * k2)
        k4 = h * func(t + h, y + k3)
        y1 = y + (k1 + 2*k2 + 2*k3 + k4)/6
        
        h_half = h * 0.5
        
        k1 = h_half * func(t, y)
        k2 = h_half * func(t + 0.5 * h_half, y + 0.5 * k1)
        k3 = h_half * func(t + 0.5 * h_half, y + 0.5 * k2)
        k4 = h_half * func(t + h_half, y + k3)
        y_half = y + (k1 + 2*k2 + 2*k3 + k4)/6

        k1 = h_half * func(t + h_half, y_half)
        k2 = h_half * func(t + 1.5 * h_half, y_half + 0.5 * k1)
        k3 = h_half * func(t + 1.5 * h_half, y_half + 0.5 * k2)
        k4 = h_half * func(t + 2 * h_half, y_half + k3)
        y2 = y_half + (k1 + 2*k2 + 2*k3 + k4)/6
        
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
    y_values = np.array(y_values, dtype=np.float32)
    
    return t_values, y_values
