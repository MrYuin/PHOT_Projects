# realistic_detector.py
import numpy as np
from scipy.signal import butter, filtfilt
from scipy.interpolate import interp1d

class DetectorOutput:
    def __init__(self, output_signal, detector_output_power, samples_per_bit, 
                 detector_SNR, quantization_error):
        self.output_signal = output_signal
        self.detector_output_power = detector_output_power
        self.samples_per_bit = samples_per_bit
        self.detector_SNR = detector_SNR
        self.quantization_error = quantization_error

def realistic_detector_sim(signal_with_noise, t, bit_duration, bit_rate, 
                           detector_threshold=0, enable_threshold=False):
    """
    真实探测器模拟（采样、量化、平滑）
    
    参数:
        signal_with_noise: 输入信号（带噪声）
        t: 时间数组
        bit_duration: 单个比特持续时间
        bit_rate: 比特率
        detector_threshold: 探测器阈值
        enable_threshold: 是否启用阈值
    
    返回:
        DetectorOutput对象，包含：
            - output_signal: 输出信号
            - detector_output_power: 探测器输出功率
            - samples_per_bit: 每比特采样点数
            - detector_SNR: 探测器等效SNR
            - quantization_error: 量化误差RMS
    """
    print('\n=== 开始真实探测器模拟 ===')
    
    # 计算时间步长和每码元采样点数
    dt = t[1] - t[0]
    samples_per_bit = int(np.round(bit_duration / dt))
    
    # 1. 光电探测（光强检测）
    detected_signal = np.abs(signal_with_noise)**2
    
    # 2. 模拟探测器带宽限制（低通滤波）
    detector_bandwidth = 0.7 * bit_rate  # 探测器带宽约为码率的%
    cutoff_normalized = detector_bandwidth / (1/dt/2)  # 归一化截止频率
    
    if cutoff_normalized < 1:
        # 设计4阶巴特沃斯低通滤波器
        b_lpf, a_lpf = butter(4, cutoff_normalized, btype='low')
        detected_signal = filtfilt(b_lpf, a_lpf, detected_signal)
        print(f'探测器带宽限制: {detector_bandwidth*1e-9:.2f} GHz (码率的 {(detector_bandwidth/bit_rate)*100:.0f}%)')
    
    # 3. ADC采样（降采样）
    adc_oversample_rate = 8  # 每个码元采样点数
    downsample_factor = max(1, int(np.round(samples_per_bit / adc_oversample_rate)))
    
    if downsample_factor > 1:
        sampled_signal = detected_signal[::downsample_factor]
        sampled_t = t[::downsample_factor]
        print(f'ADC采样率: {1/(dt*downsample_factor)*1e-9:.2f} GSa/s (每码元{adc_oversample_rate}个样点)')
    else:
        sampled_signal = detected_signal.copy()
        sampled_t = t.copy()
        print('ADC采样率: 与模拟采样率相同')
    
    # 4. 量化（模拟ADC位数）
    adc_bits = 8  # 8位ADC
    adc_levels = 2**adc_bits
    signal_min = np.min(sampled_signal)
    signal_max = np.max(sampled_signal)
    signal_range = signal_max - signal_min
    
    # 量化过程
    quantized_signal = np.round((sampled_signal - signal_min) / signal_range * (adc_levels - 1))
    quantized_signal = quantized_signal / (adc_levels - 1) * signal_range + signal_min
    
    # 添加量化噪声
    quantization_noise = sampled_signal - quantized_signal
    quantization_error_rms = np.sqrt(np.mean(quantization_noise**2))
    print(f'ADC量化位数: {adc_bits} bits, 量化误差RMS: {quantization_error_rms:.4e} W')
    
    # 5. 上采样回原始采样率（用于对比）
    if downsample_factor > 1:
        interp_func = interp1d(sampled_t, quantized_signal, kind='linear', 
                               bounds_error=False, fill_value='extrapolate')
        quantized_upsampled = interp_func(t)
    else:
        quantized_upsampled = quantized_signal.copy()
    
    # 6. 数字后处理平滑（移动平均滤波）
    smooth_window = int(np.round(samples_per_bit * 0.3))  # 平滑窗口为30%码元周期
    if smooth_window > 1:
        smooth_filter = np.ones(smooth_window) / smooth_window
        detector_output = np.convolve(quantized_upsampled, smooth_filter, mode='same')
        print(f'数字平滑窗口: {smooth_window} 样点 ({(smooth_window/samples_per_bit)*100:.1f}% 码元周期)')
    else:
        detector_output = quantized_upsampled.copy()
    
    # 7. 应用探测器最低阈值
    if enable_threshold:
        # 低于阈值的信号被截断为0
        below_threshold = detector_output < detector_threshold
        num_below_threshold = np.sum(below_threshold)
        percent_below = (num_below_threshold / len(detector_output)) * 100
        
        detector_output[below_threshold] = 0
        
        print(f'探测器阈值: {detector_threshold:.4e} W')
        print(f'  低于阈值的样点: {num_below_threshold} / {len(detector_output)} ({percent_below:.2f}%)')
    
    # 计算探测器引入的总失真
    ideal_signal = np.abs(signal_with_noise)**2
    detector_distortion = detector_output - ideal_signal
    distortion_power = np.mean(detector_distortion**2)
    signal_power_detector = np.mean(ideal_signal**2)
    detector_SNR = 10 * np.log10(signal_power_detector / distortion_power)
    
    print(f'探测器等效SNR: {detector_SNR:.2f} dB')
    print('=== 真实探测器模拟完成 ===\n')
    
    # 使用真实探测器输出重构信号
    output_signal = np.sqrt(detector_output) * np.exp(1j * np.angle(signal_with_noise))
    
    return DetectorOutput(
        output_signal=output_signal,
        detector_output_power=detector_output,
        samples_per_bit=samples_per_bit,
        detector_SNR=detector_SNR,
        quantization_error=quantization_error_rms
    )