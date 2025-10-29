# eye_diagram.py
import numpy as np
import matplotlib.pyplot as plt

def plot_eye(signal, t, bit_duration, num_bits_display=2, title='Eye Diagram'):
    """
    绘制眼图
    
    参数:
        signal: 输入信号（复数或实数）
        t: 时间数组
        bit_duration: 单个比特的持续时间
        num_bits_display: 眼图显示的比特数
        title: 图标题
    """
    # 取信号的幅度（如果是复数）
    if np.iscomplexobj(signal):
        signal = np.abs(signal)
    else:
        signal = np.real(signal)
    
    # 计算每个比特周期的采样点数
    dt = t[1] - t[0]
    samples_per_bit = int(bit_duration / dt)
    eye_window = samples_per_bit * num_bits_display
    
    # 创建图形

    plt.figure(figsize=(10, 6))
    
    # 创建眼图时间轴（归一化到比特周期）
    eye_time = np.linspace(0, num_bits_display * bit_duration, eye_window) * 1e9  # 转换为ns
    
    # 叠加多个比特周期
    num_traces = len(signal) // eye_window
    
    for i in range(num_traces):
        start_idx = i * eye_window
        end_idx = start_idx + eye_window
        
        if end_idx <= len(signal):
            trace = signal[start_idx:end_idx]
            plt.plot(eye_time, trace, 'b-', alpha=0.3, linewidth=0.5)
    
    plt.xlabel('Time (ns)')
    plt.ylabel('Amplitude')
    plt.title(title)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()