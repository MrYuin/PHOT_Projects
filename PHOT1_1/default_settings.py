import numpy as np


# 光纤传输
fibre_len = 5000
LN = 1000
dL = fibre_len/LN

alpha = 0.2*0.001;
beta2 = -1.276e-1*1e-24
beta3 = +8.119e-5*1e-36

c_light=3e8
fCenter = 193.1022e12;
fSpan = 0.5e12;
fN = 2**15;

'''Set switch of situation'''
enable_edge_delay = True
enable_rx_noise = True
enable_real_detector = True
enable_threshold = False

'''Set noise situation'''
noise_mode='ABSOLUTE' #SNR or ABSOLUTE

'''SNR case'''
SNR_dB_rx = 30

'''ABSOLUTE case'''
ABSOLUTE_rx = 1e-3

'''Detector Threshold'''
detector_threshold=0.5

'''Set of Square waves'''
bit_rate = 50e8
bit_duration = 1/bit_rate
num_bits = 128
total_time = num_bits * bit_duration

rise_time =  0.2 * bit_duration
fall_time = rise_time

binary_data = np.random.randint(0, 2, size=num_bits)



fB = fCenter
wB = 2*np.pi*fB


real_det_fc = 0.7 * bit_rate  # 截止频率 = 0.7倍比特率
