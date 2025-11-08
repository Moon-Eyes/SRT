% -----------------------------------------------------------------
% 脚本：GS_phase_retrieval.m
% 描述：使用 Gerchberg-Saxton 算法进行相位恢复，并包含
%       中间结果显示、计时和重建验证模块。
% -----------------------------------------------------------------

clc;           
clear;        
close all;  

%% 1. 初始化参数
iterative = 300;            
imagename = 'Butterfly.png';    
phaseimage = 'phase2D.png';   

iterations_to_show = [10, 20, 30, 50, 100, 200];
figure('Name', '迭代中间结果', 'NumberTitle', 'off');
plot_index = 1; 

%% 2. 准备输入数据

known_abs_spatial = imread(imagename);            
known_abs_spatial = rgb2gray(known_abs_spatial); 
known_abs_spatial = im2double(known_abs_spatial); 


[width, length] = size(known_abs_spatial);      
temp_phase = 2 * pi * rand(width, length) - pi;  
input = known_abs_spatial .* exp(1i * temp_phase); 
known_abs_fourier = abs(fft2(input));            

%% 3. GS 迭代算法
fprintf('开始 GS 算法迭代 (共 %d 次)...\n', iterative);

tic;

phase_estimate = pi * rand(width, length);

for p = 1:iterative
    
    signal_estimate_spatial = known_abs_spatial .* exp(1i * phase_estimate);
    
    temp1 = fft2(signal_estimate_spatial);
    temp_ang = angle(temp1);                           
    signal_estimate_fourier = known_abs_fourier .* exp(1i * temp_ang); 
        
    temp2 = ifft2(signal_estimate_fourier);

    phase_estimate = angle(temp2); 

   
    if ismember(p, iterations_to_show)
        temp_phase_img = phase_estimate;
        temp_phase_img(temp_phase_img < 0) = temp_phase_img(temp_phase_img < 0) + 2 * pi;
        temp_phase_img = temp_phase_img / (2 * pi);
        
        subplot(2, 3, plot_index);
        imshow(temp_phase_img);
        title(['第 ', num2str(p), ' 次迭代']);
        drawnow; 
        plot_index = plot_index + 1;
    end
end


elapsed_time = toc;
fprintf('GS 算法迭代完成。总耗时: %.2f 秒。\n', elapsed_time);

%% 4. 保存最终的相位图像

phase_estimate(phase_estimate < 0) = phase_estimate(phase_estimate < 0) + 2 * pi; 
retrieved = phase_estimate / (2 * pi); 

figure('Name', '最终结果', 'NumberTitle', 'off');
subplot(1, 2, 1);
imshow(known_abs_spatial); title('原始输入图像');
subplot(1, 2, 2);
imshow(retrieved); title('最终恢复的相位图像');
imwrite(retrieved, phaseimage);
fprintf('最终相位图像已保存为: %s\n', phaseimage);

%% 5. --- 重建验证模块 ---
% -----------------------------------------------------------------
% 该模块验证我们能否使用“保存的相位图” + “已知的傅里叶幅度”
% 来重建出“原始空间图像”。
% -----------------------------------------------------------------

fprintf('\n--- 开始重建验证模块 ---\n');

disp('请在弹出的窗口中选择您刚刚保存的相位图 (e.g., phase.png)');
[phase_file, phase_path] = uigetfile({'*.png'; '*.bmp'; '*.jpg'}, ...
                                    '请选择您刚才保存的相位图');

if isequal(phase_file, 0)
    disp('用户取消了选择。验证模块终止。');
else
    fprintf('已选择相位文件: %s\n', fullfile(phase_path, phase_file));
    
    retrieved_phase_img = imread(fullfile(phase_path, phase_file));
    retrieved_phase_img = im2double(retrieved_phase_img); 
    
    retrieved_phase_rad = retrieved_phase_img * 2 * pi;
    
    disp('... 正在使用 (已知的空间幅度 + 加载的相位) 进行衍射...');
    
    
    signal_spatial = known_abs_spatial .* exp(1i * retrieved_phase_rad);

    signal_fourier = fft2(signal_spatial);
    signal_fourier_constrained = known_abs_fourier .* exp(1i * angle(signal_fourier));

    reconstructed_spatial = ifft2(signal_fourier_constrained);
    reconstructed_image = abs(reconstructed_spatial);
    
    figure('Name', '重建验证结果', 'NumberTitle', 'off'); 
    subplot(1, 2, 1);
    imshow(known_abs_spatial);
    title('原始输入图像 (幅度)');
    
    subplot(1, 2, 2);
    imshow(reconstructed_image);
    title('从(相位图+傅里叶幅度)重建的图像');
end