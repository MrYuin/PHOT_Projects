import numpy as np
import matplotlib.pyplot as plt
import freq_set as freq
import time_set as time
import default_settings as ds
from fontTools.varLib.instancer import AxisLimits
class NoisedOutput:
    def __init__(self, t_array, SNR_var):
        self.out_t = t_array
        self.SNR = SNR_var

def calc(t_attay,tN):
    if ds.enable_rx_noise:
        signal_power_rx = np.mean(np.abs(t_attay)**2)
        if ds.noise_mode == 'ABSOLUTE':
            noise_power_rx = ds.ABSOLUTE_rx;
            SNR_rx_actual = 10 * np.log10(signal_power_rx / noise_power_rx);
        else:
            noise_power_rx = signal_power_rx / (10 ** (ds.SNR_dB_rx / 10));
            SNR_rx_actual = ds.SNR_dB_rx;

        noise_rx = np.sqrt(noise_power_rx / 2) * (np.random.randn(tN) + 1j * np.random.randn(tN))
        t_attay_with_noise = t_attay + noise_rx;
        print('We consider noise in Rx, with SNR = {:.3f} dB'.format(SNR_rx_actual) )
    else:
        t_attay_with_noise = t_attay;
        SNR_rx_actual = np.inf;
        print('We ignore noise in Rx')

    return NoisedOutput(t_attay_with_noise, SNR_rx_actual)