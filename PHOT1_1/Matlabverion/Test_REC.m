clear; format long
clight = 3e8;

%% 控制开关
enable_edge_delay = true;   % 边沿延迟开关
enable_tx_noise = false;    % 发射端噪声开关（通常关闭）
enable_rx_noise = false;        % 接收端（探测器）噪声开关
enable_realistic_detector = false;  % 真实探测器开关
enable_detector_threshold = false;   % 探测器最低阈值开关（新增）

% 噪声设置模式选择：'SNR' 或 'ABSOLUTE'
noise_mode = 'ABSOLUTE';    % 'SNR' = 相对信噪比模式, 'ABSOLUTE' = 绝对功率模式

% SNR模式参数
SNR_dB_tx = 30;             % 发射端信噪比（dB）
SNR_dB_rx = 40;             % 接收端信噪比（dB）- 探测器噪声

% 绝对噪声功率模式参数（单位：W）
noise_power_tx_absolute = 1e-3;     % 发射端绝对噪声功率 (W)
noise_power_rx_absolute = 5e-5;     % 接收端绝对噪声功率 (W) - 可设置为散粒噪声、热噪声等

% 探测器阈值参数
detector_threshold = 0.5;   % 探测器最低检测阈值 (W)，低于此值输出为0

%% 方波参数
bit_rate = 50e8; 
bit_duration = 1/bit_rate;
num_bits = 128;
total_time = num_bits * bit_duration;

% 边沿参数（上升/下降时间）
rise_time = 0.2 * bit_duration;  % 上升时间（%码元周期）
fall_time = 0.2 * bit_duration;  % 下降时间

rng(42);
binary_data = randi([0, 1], num_bits, 1);

%% 频率设置
fCenter = 193.1022e12;
lamCenter = clight/fCenter;
fSpan = 0.5e12;
fN = 2^15;
fFrom = fCenter - fSpan/2;
fTo = fCenter + fSpan/2;
df = fSpan/fN;
f = (fFrom:df:fTo-df)';
w = 2*pi*f;
lam = clight./f;

%% 时间设置
tN = fN;
tSpan = 1/df;
dt = tSpan/tN;
t = (0:tN-1)'*dt;

%% 生成带边沿延迟的方波
power_per_bit = 6;
in_t = zeros(tN, 1);

for idx = 1:tN
    t_current = t(idx);
    if t_current < total_time
        bit_index = floor(t_current / bit_duration) + 1;
        if bit_index >= 1 && bit_index <= num_bits
            bit_value = binary_data(bit_index);
            
            if enable_edge_delay
                % 计算当前位在码元内的相对时间
                t_in_bit = t_current - (bit_index - 1) * bit_duration;
                
                % 判断前一位的值
                prev_bit = 0;
                if bit_index > 1
                    prev_bit = binary_data(bit_index - 1);
                end
                
                % 上升沿：0->1
                if prev_bit == 0 && bit_value == 1 && t_in_bit < rise_time
                    in_t(idx) = (t_in_bit / rise_time) * bit_value;
                % 下降沿：1->0
                elseif prev_bit == 1 && bit_value == 0 && t_in_bit < fall_time
                    in_t(idx) = (1 - t_in_bit / fall_time) * prev_bit;
                else
                    in_t(idx) = bit_value;
                end
            else
                in_t(idx) = bit_value;
            end
        end
    end
end

% 加载波和功率
in_t = sqrt(power_per_bit) * in_t .* exp(1i*2*pi*fCenter*t);

% 添加发射端高斯白噪声（可选，通常关闭）
if enable_tx_noise
    signal_power_tx = mean(abs(in_t).^2);
    
    if strcmp(noise_mode, 'ABSOLUTE')
        % 绝对噪声功率模式
        noise_power_tx = noise_power_tx_absolute;
        SNR_tx_actual = 10*log10(signal_power_tx / noise_power_tx);
        fprintf('发射端噪声已添加 [绝对功率模式]\n');
        fprintf('  噪声功率: %.4e W, 信号功率: %.4e W\n', noise_power_tx, signal_power_tx);
        fprintf('  实际SNR: %.2f dB\n', SNR_tx_actual);
    else
        % SNR模式
        noise_power_tx = signal_power_tx / (10^(SNR_dB_tx/10));
        fprintf('发射端噪声已添加 [SNR模式]\n');
        fprintf('  目标SNR: %.1f dB, 噪声功率: %.4e W\n', SNR_dB_tx, noise_power_tx);
    end
    
    noise_tx = sqrt(noise_power_tx/2) * (randn(tN, 1) + 1i*randn(tN, 1));
    in_t = in_t + noise_tx;
end

%% 频域变换
in_f = fftshift(ifft(ifftshift(in_t)));

figure(2); subplot(211);
plot((f-fCenter)*1e-9, 10*log10(abs(in_f).^2/max(abs(in_f).^2)), 'r', 'LineWidth', 1.5);
xlabel('Frequency Offset (GHz)'); ylabel('Power (dBm)');
title('Frequency Spectrum'); xlim([-250, 250]); grid on;

subplot(212);
plot(f*1e-12, abs(in_f).^2, 'r', 'LineWidth', 1.5);
xlabel('Frequency (THz)'); ylabel('Power (a.u.)'); grid on;

%% 光纤参数
alpha = 0;%0.2*0.001;
beta2 = -1.276e-1*1e-24; 
beta3 = +8.119e-5*1e-36; 

fB = fCenter;
wB = 2*pi*fB;

%% 光纤传输
fibre_len = 50000;
LN = 1000;
dL = fibre_len/LN;

D = -alpha/2 + 1i*beta2*((w-wB).^2)/2 + 1i*beta3*((w-wB).^3)/6;

in_t0 = in_t;

fprintf('开始光纤传输模拟，光纤长度 = %.1f km\n', fibre_len/1000);
for m = 1:LN
    ou_t = ifftshift(fft(ifftshift(exp(dL*D).*(fftshift(ifft(fftshift(in_t)))))));
    in_t = ou_t;
    if any(isnan(in_t)), break; end
end
fprintf('光纤传输完成\n');

%% 接收端（探测器）添加高斯白噪声
if enable_rx_noise
    % 计算接收信号功率
    signal_power_rx = mean(abs(ou_t).^2);
    
    if strcmp(noise_mode, 'ABSOLUTE')
        % 绝对噪声功率模式
        noise_power_rx = noise_power_rx_absolute;
        SNR_rx_actual = 10*log10(signal_power_rx / noise_power_rx);
        fprintf('\n接收端噪声已添加 [绝对功率模式]\n');
        fprintf('  噪声功率: %.4e W (可理解为散粒噪声+热噪声)\n', noise_power_rx);
        fprintf('  信号功率: %.4e W\n', signal_power_rx);
        fprintf('  实际SNR: %.2f dB\n', SNR_rx_actual);
    else
        % SNR模式
        noise_power_rx = signal_power_rx / (10^(SNR_dB_rx/10));
        SNR_rx_actual = SNR_dB_rx;
        fprintf('\n接收端噪声已添加 [SNR模式]\n');
        fprintf('  目标SNR: %.1f dB\n', SNR_dB_rx);
        fprintf('  信号功率: %.4e W, 噪声功率: %.4e W\n', signal_power_rx, noise_power_rx);
    end
    
    % 生成高斯白噪声
    noise_rx = sqrt(noise_power_rx/2) * (randn(tN, 1) + 1i*randn(tN, 1));
    
    % 添加噪声到接收信号
    ou_t_with_noise = ou_t + noise_rx;
else
    ou_t_with_noise = ou_t;
    fprintf('接收端噪声未添加\n');
end

%% 真实探测器模拟（采样、量化、平滑）
if enable_realistic_detector
    fprintf('\n=== 开始真实探测器模拟 ===\n');
    
    % 计算每码元采样点数
    samples_per_bit = round(bit_duration / dt);
    
    % 1. 光电探测（光强检测）
    detected_signal = abs(ou_t_with_noise).^2;
    
    % 2. 模拟探测器带宽限制（低通滤波）
    detector_bandwidth = 0.7 * bit_rate;  % 探测器带宽约为码率的70%
    cutoff_normalized = detector_bandwidth / (1/dt/2);  % 归一化截止频率
    
    if cutoff_normalized < 1
        % 设计4阶巴特沃斯低通滤波器
        [b_lpf, a_lpf] = butter(4, cutoff_normalized);
        detected_signal = filtfilt(b_lpf, a_lpf, detected_signal);
        fprintf('探测器带宽限制: %.2f GHz (码率的 %.0f%%)\n', ...
                detector_bandwidth*1e-9, (detector_bandwidth/bit_rate)*100);
    end
    
    % 3. ADC采样（降采样）
    adc_oversample_rate = 4;  % 每个码元采样4个点
    downsample_factor = max(1, round(samples_per_bit / adc_oversample_rate));
    
    if downsample_factor > 1
        sampled_signal = detected_signal(1:downsample_factor:end);
        sampled_t = t(1:downsample_factor:end);
        fprintf('ADC采样率: %.2f GSa/s (每码元%d个样点)\n', ...
                1/(dt*downsample_factor)*1e-9, adc_oversample_rate);
    else
        sampled_signal = detected_signal;
        sampled_t = t;
        fprintf('ADC采样率: 与模拟采样率相同\n');
    end
    
    % 4. 量化（模拟ADC位数）
    adc_bits = 8;  % 8位ADC
    adc_levels = 2^adc_bits;
    signal_min = min(sampled_signal);
    signal_max = max(sampled_signal);
    signal_range = signal_max - signal_min;
    
    % 量化过程
    quantized_signal = round((sampled_signal - signal_min) / signal_range * (adc_levels - 1));
    quantized_signal = quantized_signal / (adc_levels - 1) * signal_range + signal_min;
    
    % 添加量化噪声
    quantization_noise = sampled_signal - quantized_signal;
    quantization_error_rms = sqrt(mean(quantization_noise.^2));
    fprintf('ADC量化位数: %d bits, 量化误差RMS: %.4e W\n', adc_bits, quantization_error_rms);
    
    % 5. 上采样回原始采样率（用于对比）
    if downsample_factor > 1
        quantized_upsampled = interp1(sampled_t, quantized_signal, t, 'linear', 'extrap');
    else
        quantized_upsampled = quantized_signal;
    end
    
    % 6. 数字后处理平滑（移动平均滤波）
    smooth_window = round(samples_per_bit * 0.3);  % 平滑窗口为30%码元周期
    if smooth_window > 1
        smooth_filter = ones(smooth_window, 1) / smooth_window;
        detector_output = conv(quantized_upsampled, smooth_filter, 'same');
        fprintf('数字平滑窗口: %d 样点 (%.1f%% 码元周期)\n', ...
                smooth_window, (smooth_window/samples_per_bit)*100);
    else
        detector_output = quantized_upsampled;
    end
    
    % 7. 应用探测器最低阈值（新增）
    if enable_detector_threshold
        % 低于阈值的信号被截断为0
        below_threshold = detector_output < detector_threshold;
        num_below_threshold = sum(below_threshold);
        percent_below = (num_below_threshold / length(detector_output)) * 100;
        
        detector_output(below_threshold) = 0;
        
        fprintf('探测器阈值: %.4e W\n', detector_threshold);
        fprintf('  低于阈值的样点: %d / %d (%.2f%%)\n', ...
                num_below_threshold, length(detector_output), percent_below);
    end
    
    % 计算探测器引入的总失真
    ideal_signal = abs(ou_t_with_noise).^2;
    detector_distortion = detector_output - ideal_signal;
    distortion_power = mean(detector_distortion.^2);
    signal_power_detector = mean(ideal_signal.^2);
    detector_SNR = 10*log10(signal_power_detector / distortion_power);
    
    fprintf('探测器等效SNR: %.2f dB\n', detector_SNR);
    fprintf('=== 真实探测器模拟完成 ===\n\n');
    
    % 使用真实探测器输出
    ou_t_detected = sqrt(detector_output) .* exp(1i*angle(ou_t_with_noise));
else
    detector_output = abs(ou_t_with_noise).^2;
    ou_t_detected = ou_t_with_noise;
    samples_per_bit = round(bit_duration / dt);
end

ou_f = fftshift(ifft(fftshift(ou_t_with_noise)));

%% 时域对比
figure(1);
plot(t*1e9, abs(in_t0).^2, 'b', 'LineWidth', 1.5); hold on;
plot(t*1e9, abs(ou_t_with_noise).^2, 'r', 'LineWidth', 1.5);
if enable_realistic_detector
    plot(t*1e9, detector_output, 'g--', 'LineWidth', 1.5);
end
xlabel('Time (ns)'); ylabel('Power (W)');
if enable_realistic_detector
    title('Input vs Output (with Realistic Detector)'); 
    legend('Input', 'After Fiber+Noise', 'After Realistic Detector');
elseif enable_rx_noise
    title('Input vs Output (with Detector Noise)'); 
    legend('Input', 'Output + Detector Noise');
else
    title('Input vs Output (no Detector Noise)'); 
    legend('Input', 'Output');
end
xlim([0, 5]); grid on; hold off;

%% 频域对比
figure(2);
plot((f-fCenter)*1e-9, 10*log10(abs(in_f).^2/max(abs(in_f).^2)), 'b', 'LineWidth', 1.5); hold on;
plot((f-fCenter)*1e-9, 10*log10(abs(ou_f).^2/max(abs(in_f).^2)), 'r', 'LineWidth', 1.5);
xlabel('Frequency Offset (GHz)'); ylabel('Power (dB)');
if enable_realistic_detector
    title('Spectrum Comparison (with Realistic Detector)'); 
elseif enable_rx_noise
    title('Spectrum Comparison (with Detector Noise)'); 
else
    title('Spectrum Comparison (no Detector Noise)'); 
end
legend('Input', 'Output');
xlim([-250, 250]); ylim([-80, 5]); grid on; hold off;

%% 眼图
eye_window = 4;
valid_end = min(find(t >= total_time, 1)-1, length(t));

figure(3);
if enable_realistic_detector
    ou_power = detector_output;
else
    ou_power = abs(ou_t_with_noise).^2;
end
num_traces = floor(valid_end / samples_per_bit) - eye_window + 1;

hold on;
for k = 1:num_traces
    start_idx = (k-1)*samples_per_bit + 1;
    end_idx = start_idx + samples_per_bit*eye_window - 1;
    if end_idx <= valid_end
        trace_time = ((0:samples_per_bit*eye_window-1) - samples_per_bit) * dt * 1e9;
        plot(trace_time, ou_power(start_idx:end_idx), 'k', 'LineWidth', 0.5);
    end
end
xlabel('Time (ns)'); ylabel('Power (W)');
if enable_realistic_detector
    title('Eye Diagram (Realistic Detector)'); 
elseif enable_rx_noise 
    title(sprintf('Eye Diagram (Detector SNR = %.1f dB)', SNR_dB_rx)); 
else
    title('Eye Diagram (No Detector Noise)');
end
xlim([-(bit_duration*1e9), (2*bit_duration*1e9)]); 
grid on;
yrange = ylim;
plot([0, 0], yrange, 'r--', 'LineWidth', 1.5);
hold off;

%% 探测器性能详细分析（新增）
if enable_realistic_detector
    figure(5);
    
    % 子图1：探测器各阶段信号对比
    subplot(221);
    plot(t*1e9, abs(ou_t_with_noise).^2, 'b', 'LineWidth', 1); hold on;
    plot(t*1e9, detected_signal, 'r', 'LineWidth', 1);
    plot(t*1e9, detector_output, 'g', 'LineWidth', 1.5);
    if enable_detector_threshold
        plot([0, max(t)*1e9], [detector_threshold, detector_threshold], ...
             'k--', 'LineWidth', 2);
        legend('Ideal (PD input)', 'After BW limit', 'After ADC+Filter+Threshold', 'Detector Threshold');
    else
        legend('Ideal (PD input)', 'After BW limit', 'After ADC+Filter');
    end
    xlabel('Time (ns)'); ylabel('Power (W)');
    title('Detector Signal Processing Stages');
    xlim([0, 5]); grid on; hold off;
    
    % 子图2：阈值效应分析
    subplot(222);
    if enable_detector_threshold
        % 计算阈值前后的信号
        detector_before_threshold = detector_output;
        if smooth_window > 1
            detector_before_threshold = conv(quantized_upsampled, smooth_filter, 'same');
        else
            detector_before_threshold = quantized_upsampled;
        end
        
        % 显示被阈值截断的部分
        plot(t*1e9, detector_before_threshold, 'b', 'LineWidth', 1); hold on;
        plot(t*1e9, detector_output, 'r', 'LineWidth', 1.5);
        plot([0, max(t)*1e9], [detector_threshold, detector_threshold], ...
             'k--', 'LineWidth', 2);
        xlabel('Time (ns)'); ylabel('Power (W)');
        title('Threshold Effect');
        legend('Before Threshold', 'After Threshold', sprintf('Threshold (%.2e W)', detector_threshold));
        xlim([0, 5]); grid on; hold off;
    else
        % 显示量化误差
        if downsample_factor > 1
            plot(sampled_t*1e9, quantization_noise, 'k', 'LineWidth', 0.5);
            xlabel('Time (ns)'); ylabel('Quantization Error (W)');
            title(sprintf('ADC Quantization Noise (%d bits)', adc_bits));
            xlim([0, 5]); grid on;
        else
            text(0.5, 0.5, 'No downsampling performed', ...
                 'HorizontalAlignment', 'center', 'Units', 'normalized');
            title('Quantization Analysis');
        end
    end
    
    % 子图3：信号功率谱密度对比
    subplot(223);
    ideal_psd = pwelch(abs(ou_t_with_noise).^2, [], [], [], 1/dt);
    detector_psd = pwelch(detector_output, [], [], [], 1/dt);
    freq_psd = linspace(0, 1/dt/2, length(ideal_psd));
    plot(freq_psd*1e-9, 10*log10(ideal_psd), 'b', 'LineWidth', 1); hold on;
    plot(freq_psd*1e-9, 10*log10(detector_psd), 'r', 'LineWidth', 1.5);
    xlabel('Frequency (GHz)'); ylabel('PSD (dB/Hz)');
    title('Power Spectral Density');
    legend('Ideal', 'Realistic Detector');
    xlim([0, bit_rate*1e-9*2]); grid on; hold off;
    
    % 子图4：信号幅度分布直方图
    subplot(224);
    histogram(abs(ou_t_with_noise).^2, 50, 'Normalization', 'probability', ...
              'FaceColor', 'b', 'FaceAlpha', 0.5, 'EdgeColor', 'none'); hold on;
    histogram(detector_output, 50, 'Normalization', 'probability', ...
              'FaceColor', 'r', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    if enable_detector_threshold
        plot([detector_threshold, detector_threshold], ylim, 'k--', 'LineWidth', 2);
        legend('Ideal Signal', 'Detector Output', 'Threshold');
    else
        legend('Ideal Signal', 'Detector Output');
    end
    xlabel('Power (W)'); ylabel('Probability');
    title('Signal Amplitude Distribution');
    grid on; hold off;
    
    sgtitle('Realistic Detector Performance Analysis');
end

%% 信号质量统计
if enable_rx_noise || enable_realistic_detector
    figure(4);
    subplot(211);
    plot(t*1e9, abs(ou_t).^2, 'b', 'LineWidth', 1); hold on;
    plot(t*1e9, abs(ou_t_with_noise).^2, 'r', 'LineWidth', 1);
    if enable_realistic_detector
        plot(t*1e9, detector_output, 'g--', 'LineWidth', 1.5);
        legend('Ideal (no noise)', 'With Rx Noise', 'With Realistic Detector');
    else
        legend('Before Detector', 'After Detector');
    end
    xlabel('Time (ns)'); ylabel('Power (W)');
    title('Signal Comparison: Processing Stages');
    xlim([0, 5]); grid on; hold off;
    
    subplot(212);
    if enable_realistic_detector
        total_degradation = detector_output - abs(ou_t).^2;
        plot(t*1e9, abs(total_degradation), 'k', 'LineWidth', 0.5);
        title('Total Signal Degradation (Noise + Detector)');
    else
        noise_only = abs(ou_t_with_noise - ou_t).^2;
        plot(t*1e9, noise_only, 'k', 'LineWidth', 0.5);
        title('Detector Noise Power');
    end
    xlabel('Time (ns)'); ylabel('Degradation Power (W)');
    xlim([0, 5]); grid on;
end

%% 输出总结
fprintf('\n=== 系统性能总结 ===\n');
fprintf('码率: %.2f Gb/s\n', bit_rate*1e-9);
fprintf('光纤长度: %.1f km\n', fibre_len/1000);
fprintf('噪声模式: %s\n', noise_mode);
if enable_rx_noise
    if strcmp(noise_mode, 'ABSOLUTE')
        fprintf('接收端噪声功率: %.4e W (实际SNR: %.2f dB)\n', ...
                noise_power_rx, SNR_rx_actual);
    else
        fprintf('接收端SNR: %.1f dB\n', SNR_dB_rx);
    end
end
if enable_realistic_detector
    fprintf('探测器等效SNR: %.2f dB\n', detector_SNR);
    fprintf('ADC: %d bits, 采样率 %.2f GSa/s\n', adc_bits, 1/(dt*downsample_factor)*1e-9);
    if enable_detector_threshold
        fprintf('探测器阈值: %.4e W\n', detector_threshold);
    end
end
fprintf('==================\n');