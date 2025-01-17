import numpy as np

from aux.var import A, B, XI
from aux.func import func
from aux.init_h import init_h

S = 2
def rk2wc_cpu(y0, t0, t_end, tol, flag=0):

    a21 = XI
    b2 = 1 / (2 * XI)
    b1 = 1 - b2
    
    t = t0
    y = y0
    t_values = [t]
    y_values = [y]

    h = init_h(y0, t0, t_end, tol, S)
    
    if flag in (1, 4):
        ta, ha, ea = [], [], []

    cnt = 0
    while t < t_end:
        if t + h > t_end:
            h = t_end - t

        k1 = h * func(t, y)
        k2 = h * func(t + a21 * h, y + a21 * k1)
        y1 = y + b1 * k1 + b2 * k2
        
        h_half = h * 0.5
        k1 = h_half * func(t, y)
        k2 = h_half * func(t + a21 * h_half, y + a21 * k1)
        y_half = y + b1 * k1 + b2 * k2

        t += h_half

        k1 = h_half * func(t, y_half)
        k2 = h_half * func(t + a21 * h_half, y_half + a21 * k1)
        y2 = y_half + b1 * k1 + b2 * k2
        
        cnt += 5
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

        if flag in (1, 4):
            ta.append(t)
            ha.append(h)
            ea.append(error)
            
    t_values = np.array(t_values)
    y_values = np.array(y_values, dtype=np.float64)

    if flag == 1:
        return ta, ha
    elif flag == 2:
        return cnt
    elif flag == 4:
        return ta, ha, ea
    return t_values, y_values
