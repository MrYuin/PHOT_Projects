import numpy as np

def set(t_N,t_Span):
    dt=t_Span/t_N
    t = np.arange(0,t_N) * dt
    return t
