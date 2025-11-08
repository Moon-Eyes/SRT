% -----------------------------------------------------------------
% GS_ASM.m
% 结合角谱法(ASM)的GS算法，计算菲涅尔衍射全息图
% -----------------------------------------------------------------

clc;   
close all;

%% --- 1. 定义物理参数 (参考 circle_diffraction.m) ---

lambda = 532e-9; 
z = 10e-3;      
N = 1024;         

spot_diameter_physical = 0.5e-3; 
spot_to_image_ratio = 0.5;
field_size_physical = spot_diameter_physical / spot_to_image_ratio; 
dx = field_size_physical / N; 

fprintf('--------- 仿真物理参数 ---------\n');
fprintf('波长 (lambda):       %.3f nm\n', lambda * 1e9);
fprintf('衍射距离 (z):        %.3f mm\n', z * 1e3);
fprintf('仿真区域物理边长 (L): %.3f mm\n', field_size_physical * 1e3);
fprintf('单个像素物理尺寸 (dx): %.3f µm \n', dx * 1e6);
fprintf('总像素数 (N x N):    %d x %d\n', N, N);
fprintf('----------------------------------\n\n');

%% --- 2. 计算角谱衍射的传递函数 (H) ---

% 空间频率网格
df = 1 / field_size_physical; 
fx = ((-N/2 : N/2-1) * df);  
fy = ((-N/2 : N/2-1) * df);  
[Fx, Fy] = meshgrid(fx, fy); %

% 传递函数 H (用于正向传播 +z)
k_term_squared = (1/lambda)^2 - Fx.^2 - Fy.^2; %
kz = sqrt(k_term_squared);
H = exp(1j * 2 * pi * z * kz); %

% 反向传递函数 H_inv (用于反向传播 -z)
% 也就是 H 的复共轭
H_inv = conj(H);


%% --- 3. 载入“目标”：衍射平面的“已知幅度” ---

iterative = 300;                % 迭代次数
imagename = 'Butterfly.png';    % 目标图像
phaseimage = 'phase hologram ASM.png'; % 要保存的相位全息图

target_abs_z = imread(imagename);
target_abs_z = rgb2gray(target_abs_z);
target_abs_z = im2double(target_abs_z);
target_abs_z = imresize(target_abs_z, [N, N]); % 调整大小以匹配仿真网格

% 注意：我们不再需要 fftshift 目标图像，
% 因为ASM是在真实空间坐标系中工作的。

%% --- 4. “物平面”的“已知幅度”（约束） ---

% 振幅恒定为 1 (模拟平面波入射)
known_abs_spatial = ones(N, N);

%% --- 5. 初始化：物平面的随机相位 ---

% 初始相位是 [0, 2*pi] 之间的随机值
phase_estimate_spatial = 2 * pi * rand(N, N);

%% --- 6. GS 迭代循环 (使用角谱法 ASM) ---

fprintf('开始GS迭代 (使用ASM)...\n');
for p = 1:iterative
    
    % --- 步骤 1: 构造“物平面”复数场 (z=0) ---
    % (施加“物平面”幅度约束: 振幅为1)
    signal_spatial_U0 = known_abs_spatial .* exp(1i * phase_estimate_spatial);
    
    % --- 步骤 2: 传播到“衍射平面” (ASM 正向传播 +z) ---
    U0_freq = fftshift(fft2(signal_spatial_U0)); %
    Uz_freq = U0_freq .* H; %
    signal_z_Uz = ifft2(ifftshift(Uz_freq)); %
    
    % --- 步骤 3: 施加“衍射平面”幅度约束 (z=z) ---
    % (用目标图像 'Butterfly.png' 的幅度替换，保留计算出的相位)
    temp_ang_z = angle(signal_z_Uz);
    signal_z_constrained = target_abs_z .* exp(1i * temp_ang_z);
    
    % --- 步骤 4: 传播回“物平面” (ASM 反向传播 -z) ---
    Uz_constrained_freq = fftshift(fft2(signal_z_constrained));
    U0_new_freq = Uz_constrained_freq .* H_inv; % 使用反向传递函数
    signal_spatial_new = ifft2(ifftshift(U0_new_freq));
    
    % --- 步骤 5: 提取“物平面”的新相位 ---
    % (振幅被丢弃，将在下一次循环的步骤1中被重置为1)
    phase_estimate_spatial = angle(signal_spatial_new);
    
    if mod(p, 50) == 0
        fprintf('迭代 %d / %d\n', p, iterative);
    end
end
fprintf('迭代完成。\n');

%% --- 7. 结果处理和显示 ---

% 'phase_estimate_spatial' 就是我们最终需要的相位全息图
% 将其从 [-pi, +pi] 映射到 [0, 2pi]
phase_estimate_spatial(phase_estimate_spatial < 0) = phase_estimate_spatial(phase_estimate_spatial < 0) + 2 * pi;
% 归一化到 [0, 1] 以便保存为图像
retrieved_phase_hologram = phase_estimate_spatial / (2 * pi);

figure('Name', 'GS-ASM 结果')
subplot(1, 2, 1);
imshow(retrieved_phase_hologram);
title('计算出的相位全息图 (物平面, z=0)')
imwrite(retrieved_phase_hologram, phaseimage)

%% --- 8. (可选) 验证重建效果 ---

% 用平面波(振幅=1)照射我们计算出的全息图 (在 z=0)
final_object_field = 1 .* exp(1i * phase_estimate_spatial);

% 传播到衍射平面 (z=z)
final_U0_freq = fftshift(fft2(final_object_field));
final_Uz_freq = final_U0_freq .* H;
reconstruction_at_z = ifft2(ifftshift(final_Uz_freq));

% 查看重建图像的 *强度* (幅度的平方)
reconstruction_intensity = abs(reconstruction_at_z).^2;

% 归一化并显示
reconstruction_display = reconstruction_intensity / max(reconstruction_intensity(:));

subplot(1, 2, 2);
x_coords_mm = ((-N/2 : N/2-1) * dx) * 1e3; % 物理坐标 (mm)
y_coords_mm = ((-N/2 : N/2-1) * dx) * 1e3; % 物理坐标 (mm)
imagesc(x_coords_mm, y_coords_mm, reconstruction_display);
axis square;
colormap(gray);
title(['模拟重建图像 (衍射平面, z=', num2str(z*1e3), ' mm)']);
xlabel('X 坐标 (mm)');
ylabel('Y 坐标 (mm)');