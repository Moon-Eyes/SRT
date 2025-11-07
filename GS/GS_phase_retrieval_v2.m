% -----------------------------------------------------------------
% 脚本：GS_phase_retrieval_v2.m
% 描述：使用 Gerchberg-Saxton 算法进行相位恢复，并包含
%       中间结果显示、计时和重建验证模块。
% -----------------------------------------------------------------

clc;           
clear;        
close all;  

%% 1. 初始化参数
iterative = 300;            % 设迭代次数为300次
imagename = 'Butterfly.png';    % 想要提取相位的图像名称 (已为您修复)
phaseimage = 'phase2.png';   % 要保存的相位图像名称

% --- 新增功能：设置要显示的中间迭代次数 ---
iterations_to_show = [1, 2, 3, 5, 10, 20];
figure('Name', '迭代中间结果', 'NumberTitle', 'off'); % 为中间结果创建一个新窗口
plot_index = 1; % 子图索引

%% 2. 准备输入数据
% 空域输入图像的幅度（是已知的，也就是清晰的图像，它的灰度就是幅值）
known_abs_spatial = imread(imagename);            % 作为输入图像的幅度，是已知的
known_abs_spatial = rgb2gray(known_abs_spatial); % 转换为灰度图像
known_abs_spatial = im2double(known_abs_spatial); % 将图像灰度映射到0～1

% 创建一个“目标”傅里叶幅度
% (在实际应用中，这可能是已知的衍射图样。这里我们用原始图像加一个随机相位来模拟生成它)
[width, length] = size(known_abs_spatial);       % 获取图像大小
temp_phase = 2 * pi * rand(width, length) - pi;  % [-pi, +pi] 范围内的随机相位
input = known_abs_spatial .* exp(1i * temp_phase); % 最终输入图像:幅度*e^(i*相位角度)
known_abs_fourier = abs(fft2(input));            % 傅氏变换后的幅度 (我们的“已知”约束)

%% 3. GS 迭代算法
fprintf('开始 GS 算法迭代 (共 %d 次)...\n', iterative);

% --- 新增功能：开始计时 ---
tic;

% 初始相位估计：使用随机相位
phase_estimate = pi * rand(width, length); % 像素值在[0,pi]范围内随机生成

% --- 开始迭代 ---
for p = 1:iterative
    % Step 1: 构造空间信号 (施加空间幅度约束)
    signal_estimate_spatial = known_abs_spatial .* exp(1i * phase_estimate);
    
    % Step 2: 傅立叶变换到频域，并施加傅里叶幅度约束
    temp1 = fft2(signal_estimate_spatial);
    temp_ang = angle(temp1);                           % 保留计算出的傅里叶相位
    signal_estimate_fourier = known_abs_fourier .* exp(1i * temp_ang); % 替换为已知的傅里叶幅度
    
    % Step 3: 傅立叶反变换回空域
    temp2 = ifft2(signal_estimate_fourier);
    
    % Step 4: 提取相位，作为下一次迭代的估计 (施加空间幅度约束在下一次循环的Step 1)
    phase_estimate = angle(temp2); % 范围是 [-pi, pi]

    % --- 新增功能：显示特定的中间结果 ---
    if ismember(p, iterations_to_show)
        % 将相位从 [-pi, pi] 转换到 [0, 1] 以便显示
        temp_phase_img = phase_estimate;
        temp_phase_img(temp_phase_img < 0) = temp_phase_img(temp_phase_img < 0) + 2 * pi;
        temp_phase_img = temp_phase_img / (2 * pi);
        
        subplot(2, 3, plot_index);
        imshow(temp_phase_img);
        title(['第 ', num2str(p), ' 次迭代']);
        drawnow; % 立即刷新图像窗口
        plot_index = plot_index + 1;
    end
end

% --- 新增功能：停止计时并报告时间 ---
elapsed_time = toc;
fprintf('GS 算法迭代完成。总耗时: %.2f 秒。\n', elapsed_time);

%% 4. 保存最终的相位图像
% 把 estimate_phase 从 [-pi,+pi]，映射到 [0, 2pi]
phase_estimate(phase_estimate < 0) = phase_estimate(phase_estimate < 0) + 2 * pi; 
retrieved = phase_estimate / (2 * pi); % 再映射到 [0, 1]

figure('Name', '最终结果', 'NumberTitle', 'off');
subplot(1, 2, 1);
imshow(known_abs_spatial); title('原始输入图像');
subplot(1, 2, 2);
imshow(retrieved); title('最终恢复的相位图像');
imwrite(retrieved, phaseimage);
fprintf('最终相位图像已保存为: %s\n', phaseimage);




%% 5. --- 新增功能：重建验证模块 ---
% -----------------------------------------------------------------
% 该模块验证我们能否使用“保存的相位图” + “已知的傅里叶幅度”
% 来重建出“原始空间图像”。
% -----------------------------------------------------------------

fprintf('\n--- 开始重建验证模块 ---\n');
% 提示用户选择刚刚保存的相位文件
disp('请在弹出的窗口中选择您刚刚保存的相位图 (e.g., phase.png)');
[phase_file, phase_path] = uigetfile({'*.png'; '*.bmp'; '*.jpg'}, ...
                                    '请选择您刚才保存的相位图');

if isequal(phase_file, 0)
    disp('用户取消了选择。验证模块终止。');
else
    fprintf('已选择相位文件: %s\n', fullfile(phase_path, phase_file));
    
    % 1. 加载并处理保存的相位图
    retrieved_phase_img = imread(fullfile(phase_path, phase_file));
    retrieved_phase_img = im2double(retrieved_phase_img); % 转换为 [0, 1]
    
    % 将 [0, 1] 映射回 [0, 2*pi] 的相位弧度
    % 注意: exp(1i*phase) 是 2*pi 周期的, 所以 [0, 2*pi] 和 [-pi, pi] 等效
    retrieved_phase_rad = retrieved_phase_img * 2 * pi;
    
    % 2. 执行重建 (模拟衍射过程)
    % 我们需要工作区中的原始空间幅度 (known_abs_spatial)
    % 和 目标傅里叶幅度 (known_abs_fourier)
    
    disp('... 正在使用 (已知的空间幅度 + 加载的相位) 进行衍射...');
    
    % Step A: 构造空间复振幅 (使用我们的相位)
    % 这是模拟一个空间光调制器 (SLM) 加载了我们的相位，并被原始图像形状的光束照射
    signal_spatial = known_abs_spatial .* exp(1i * retrieved_phase_rad);
    
    % Step B: 衍射到傅里叶平面 (FFT)
    signal_fourier = fft2(signal_spatial);
    
    % Step C: 施加傅里叶平面的幅度约束 (用已知的傅里叶幅度替换计算出的幅度)
    % 这模拟了光波通过一个只允许特定衍射图样通过的“滤波器”
    signal_fourier_constrained = known_abs_fourier .* exp(1i * angle(signal_fourier));
    
    % Step D: 衍射回空间平面 (IFFT)
    reconstructed_spatial = ifft2(signal_fourier_constrained);
    
    % 3. 提取重建后的幅度（这就是重建的图像）
    reconstructed_image = abs(reconstructed_spatial);
    
    % 4. 显示对比结果
    figure('Name', '重建验证结果', 'NumberTitle', 'off');
    subplot(1, 2, 1);
    imshow(known_abs_spatial);
    title('原始输入图像 (幅度)');
    
    subplot(1, 2, 2);
    imshow(reconstructed_image);
    title('从(相位图+傅里叶幅度)重建的图像');
    
    fprintf('--- 验证完成 --- \n');
end