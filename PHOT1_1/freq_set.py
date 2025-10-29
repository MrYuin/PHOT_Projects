import numpy as np

def freq_set(c_light,f_center,f_span,f_N):
    lam_center=c_light/f_center
    f_From = f_center - f_span/2
    f_To = f_center + f_span/2
    df = f_span/f_N
    f = np.arange(f_From, f_To, df)
    return f
