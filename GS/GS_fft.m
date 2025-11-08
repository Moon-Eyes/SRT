iterative=300;
imagename='Butterfly.png';    % 目标图像 (在傅里叶平面)
phaseimage='phase_hologram.png'; % 要保存的相位全息图 (在物平面)

% ------------------------------------------------
% 1. 载入“目标”：傅里叶平面的“已知幅度”
% ------------------------------------------------
target_abs_fourier = imread(imagename);
target_abs_fourier = rgb2gray(target_abs_fourier);
target_abs_fourier = im2double(target_abs_fourier);

% 傅里叶变换的零频在角落(1,1)，而我们的目标图像中心在中间
% 我们需要用 fftshift 将图像中心移到(1,1)，使其与fft2的输出对齐
target_abs_fourier = fftshift(target_abs_fourier);

[width, length] = size(target_abs_fourier);

% ------------------------------------------------
% 2. “物平面”的“已知幅度”（约束）
% ------------------------------------------------
% 振幅恒定为 1 (模拟平面波入射)
known_abs_spatial = ones(width, length);

% ------------------------------------------------
% 3. 初始化：物平面的随机相位
% ------------------------------------------------
% 初始相位是 [0, 2*pi] 之间的随机值
phase_estimate_spatial = 2 * pi * rand(width, length);

% ------------------------------------------------
% 4. GS 迭代循环 (按照您的思路)
% ------------------------------------------------
for p = 1:iterative
    % 步骤 1: 构造“物平面”复数场
    % (施加“物平面”幅度约束: 振幅为1)
    signal_spatial = known_abs_spatial .* exp(1i * phase_estimate_spatial);
    
    % 步骤 2: 传播到“傅里叶平面” (FFT)
    signal_fourier = fft2(signal_spatial);
    
    % 步骤 3: 施加“傅里叶平面”幅度约束
    % (用目标图像 'Butterfly.png' 的幅度替换，保留计算出的相位)
    temp_ang_fourier = angle(signal_fourier);
    signal_fourier_constrained = target_abs_fourier .* exp(1i * temp_ang_fourier);
    
    % 步骤 4: 传播回“物平面” (IFFT)
    signal_spatial_new = ifft2(signal_fourier_constrained);
    
    % 步骤 5: 提取“物平面”的新相位
    % (振幅被丢弃，将在下一次循环的步骤1中被重置为1)
    phase_estimate_spatial = angle(signal_spatial_new);
end
% ------------------------------------------------
% 5. 结果处理和显示
% ------------------------------------------------
% 'phase_estimate_spatial' 就是我们最终需要的相位全息图
% 将其从 [-pi, +pi] 映射到 [0, 2pi]
phase_estimate_spatial(phase_estimate_spatial < 0) = phase_estimate_spatial(phase_estimate_spatial < 0) + 2 * pi;
% 归一化到 [0, 1] 以便保存为图像
retrieved_phase_hologram = phase_estimate_spatial / (2 * pi);

figure(1)
imshow(retrieved_phase_hologram);
title('计算出的相位全息图（物平面）')
imwrite(retrieved_phase_hologram, phaseimage)


% ------------------------------------------------
% 6. (可选) 验证重建效果
% ------------------------------------------------
% 让我们看看我们计算出的全息图是否真的能重建出蝴蝶
% 用平面波(振幅=1)照射我们计算出的全息图
final_object_field = 1 .* exp(1i * phase_estimate_spatial);
% 传播到傅里叶平面
reconstruction_fourier = fft2(final_object_field);
% 查看重建图像的 *强度* (幅度的平方)
reconstruction_intensity = abs(reconstruction_fourier).^2;

% 由于我们之前对目标做了 fftshift,
% 现在需要 ifftshift 才能将重建图像的中心移回视觉中心
reconstruction_display = ifftshift(reconstruction_intensity);

% 归一化并用对数尺度显示，以便看清细节 (原始强度对比度可能很高)
reconstruction_display_log = log(1 + reconstruction_display);
reconstruction_display_log = mat2gray(reconstruction_display_log);

figure(2)
imshow(reconstruction_display_log);
title('模拟重建的傅里叶平面图像')