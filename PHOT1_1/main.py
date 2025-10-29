import numpy as np
import matplotlib.pyplot as plt
import freq_set as freq
import time_set as time
import default_settings as ds
import noise_calc as noise
import eye_diagram as eye
import scipy as scipy
import realistic_detector as rd
from scipy.signal import butter, lfilter
from fontTools.varLib.instancer import AxisLimits
from numpy.f2py.crackfortran import endifs

# (c_light, f_center, f_span, f_N)
f = freq.freq_set(ds.c_light, ds.fCenter, ds.fSpan, ds.fN)

'''Creat the square wave in time domain'''
tN=ds.fN
t_Span=ds.fN/ds.fSpan
# set(t_N, t_Span)
t = time.set(tN,t_Span)


# 生成带边沿延迟的方波
power_per_bit = 6
in_t = np.zeros(tN)

# 找到有效时间范围内的索引
valid_mask = t < ds.total_time
valid_indices = np.where(valid_mask)[0]


for idx in valid_indices:
    t_current = t[idx]
    bit_index = int(np.floor(t_current / ds.bit_duration))

    if 0 <= bit_index < ds.num_bits:
        bit_value = ds.binary_data[bit_index]

        if ds.enable_edge_delay:
            t_in_bit = t_current - bit_index * ds.bit_duration
            prev_bit = 0 if bit_index == 0 else ds.binary_data[bit_index - 1]

            # rise：0->1
            if prev_bit == 0 and bit_value == 1 and t_in_bit < ds.rise_time:
                in_t[idx] = (t_in_bit / ds.rise_time) * bit_value
            # fall：1->0
            elif prev_bit == 1 and bit_value == 0 and t_in_bit < ds.fall_time:
                in_t[idx] = (1 - t_in_bit / ds.fall_time) * prev_bit
            else:
                in_t[idx] = bit_value
        else:
            in_t[idx] = bit_value

in_t = np.sqrt(power_per_bit) * in_t * np.exp(1j*2*np.pi*ds.fCenter*t);
'''Creat square wave in time domain'''


in_f = np.fft.fftshift(np.fft.ifft(np.fft.ifftshift(in_t)));

w=2*np.pi*f

D = -ds.alpha/2 + 1j*ds.beta2*((w-ds.wB)**2)/2 + 1j*ds.beta3*((w-ds.wB)**3)/6
in_t0=in_t
in_f0=in_f

print(f'Start the propagation in fiber of = {ds.fibre_len/1000:.1f} km')
for m in range(ds.LN):
    out_t = np.fft.ifftshift(np.fft.fft(np.fft.ifftshift(np.exp(ds.dL*D) * (np.fft.fftshift(np.fft.ifft(np.fft.fftshift(in_t)))))))
    in_t = out_t
    if np.any(np.isnan(in_t)):
        break
print('Finish')




out_f = np.fft.fftshift(np.fft.ifft(np.fft.ifftshift(out_t)));


eye.plot_eye(out_t, t, ds.bit_duration, num_bits_display=3, 
             title=f'Output Eye Diagram Ideal Detector (No Noise)')


out_t_noise = noise.calc(out_t,tN);
noised_out_t = out_t_noise.out_t
SNR = out_t_noise.SNR

if ds.enable_real_detector:
    detector_result = rd.realistic_detector_sim(
        noised_out_t, t, ds.bit_duration, ds.bit_rate,
        ds.detector_threshold, ds.enable_threshold
    )
    final_signal = detector_result.output_signal
    title_info = f'Real Detector, Detector SNR={detector_result.detector_SNR:.2f} dB'
else:
    final_signal = noised_out_t
    title_info = f'SNR = {SNR:.2f} dB' if ds.enable_rx_noise else 'No Noise'

# 然后绘制
eye.plot_eye(final_signal, t, ds.bit_duration, num_bits_display=5,
             title=f'Output Eye Diagram ({title_info})')


plt.plot(f, 10*np.log10(np.abs(in_f)), label='Input Signal',color='blue')
plt.plot(f, 10*np.log10(np.abs(out_f)), label='Output Signal',color='red')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Magnitude (dB)')
plt.title('Frequency Domain - Magnitude Spectrum')
plt.grid(True)
plt.legend()
plt.show()