% ------------------------------------------------
% 1. 让用户选择一个相位全息图文件
% ------------------------------------------------
[fileName, filePath] = uigetfile({'*.png';'*.bmp';'*.jpg';'*.tif'}, '请选择一个归一化的相位图像');

% 检查用户是否取消了选择
if isequal(fileName, 0)
    disp('用户取消了操作');
    return;
end

% 拼接完整的文件路径
fullFileName = fullfile(filePath, fileName);

% ------------------------------------------------
% 2. 加载图像并将其转换回相位
% ------------------------------------------------
% 读取保存的相位图像
% 这是一个归一化的 [0, 1] 范围内的灰度图
normalized_phase_image = imread(fullFileName);

% 如果图像是彩色的（例如某些png保存格式），先转为灰度
if size(normalized_phase_image, 3) > 1
    normalized_phase_image = rgb2gray(normalized_phase_image);
end

% 将其转换为 double 类型
normalized_phase_image = im2double(normalized_phase_image);

% “反归一化”
% 之前保存时是： (相位值 [0, 2*pi]) / (2*pi) = 图像 [0, 1]
% 现在反过来：   (图像 [0, 1]) * (2*pi) = 相位值 [0, 2*pi]
phase_map_rad = normalized_phase_image * 2 * pi;

% ------------------------------------------------
% 3. 模拟重建过程
% ------------------------------------------------

% 模拟一个振幅恒为 1 的平面波照射到这个相位板上
% 'final_object_field' 是物平面的复数场
final_object_field = 1 .* exp(1i * phase_map_rad);

% 传播到傅里叶平面 (即进行傅里叶变换)
reconstruction_fourier = fft2(final_object_field);

% 计算傅里叶平面的 *强度* (Intensity)
% 强度是幅度的平方
reconstruction_intensity = abs(reconstruction_fourier).^2;

% ------------------------------------------------
% 4. 显示重建结果
% ------------------------------------------------

% 使用 ifftshift 将傅里叶变换的零频（中心）从角落移回图像中心，以便查看
reconstruction_display = ifftshift(reconstruction_intensity);

% 傅里叶平面的强度通常对比度非常高（中心点很亮）
% 使用对数尺度 (log scale) 可以更好地观察到暗处的细节
reconstruction_display_log = log(1 + reconstruction_display);

% 将对数尺度后的图像归一化到 [0, 1] 范围以便显示
reconstruction_display_log = mat2gray(reconstruction_display_log);

% 显示最终的重建图像
figure
imshow(reconstruction_display_log);
title(['从 "', fileName, '" 重建的图像']);