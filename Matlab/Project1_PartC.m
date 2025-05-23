%% ESE5310 Project 1 - Part C:  MULTI-CHANNEL FIR FILTER-BANK IN TWO DIMENSIONS
% Student Name: [Hesham Nabil Maher]
clear all;clear all;clc;
%% Part C(a): 1D 2-Channel Perfect Reconstruction Filter Bank
load('nspeech2.mat');
x = nspeech2(:)';  

% Design a 2-channel PR-FB 
N = 9;  
[h0, h1, g0, g1] = firpr2chfb(N, 0.45);  % Band-edge at ~0.45π

% Analysis (decomposition)
x0 = filter(h0, 1, x);  % Lowpass
x1 = filter(h1, 1, x);  % Highpass

% Downsample (decimation)
x0_d = x0(1:2:end);
x1_d = x1(1:2:end);

% Upsample (interpolation) & Synthesis
x0_u = zeros(1, 2*length(x0_d)); x0_u(1:2:end) = x0_d;
x1_u = zeros(1, 2*length(x1_d)); x1_u(1:2:end) = x1_d;

y0 = filter(g0, 1, x0_u);  % Lowpass reconstruction
y1 = filter(g1, 1, x1_u);  % Highpass reconstruction

% Combined reconstruction
y = y0 + y1;

% filter delay
delay = N;  % Group delay of FIR filters
y = y(delay+1:end);
x_aligned = x(1:length(y));  % Trim original to match y

% Compute reconstruction error
error = x_aligned - y;
SNR = 20*log10(norm(x_aligned)/norm(error));
disp(['Reconstruction SNR: ', num2str(SNR), ' dB']);



figure;
% 1. Time Domain: Original vs Reconstructed (Linear)
subplot(3,1,1);
plot(x_aligned, 'b'); hold on; 
plot(y, 'r--'); 
legend('Original', 'Reconstructed');
title('Time Domain: Original vs Reconstructed');
xlabel('Samples'); ylabel('Amplitude');
grid on;

% 2. Error Signal (Linear)
subplot(3,1,2);
plot(error, 'g');
title('Reconstruction Error (Linear Scale)');
xlabel('Samples'); ylabel('Error');
grid on;

% 3. Frequency Domain (dB Scale)
NFFT = 1024;
subplot(3,1,3);
[P_orig, f] = pwelch(x_aligned, hamming(NFFT), NFFT/2, NFFT, 1);
[P_recon, ~] = pwelch(y, hamming(NFFT), NFFT/2, NFFT, 1);
plot(f, 10*log10(P_orig), 'b'); hold on;
plot(f, 10*log10(P_recon), 'r--');
title('Frequency Domain (dB Scale)');
xlabel('Normalized Frequency (×π rad/sample)'); 
ylabel('Magnitude (dB)');
legend('Original', 'Reconstructed');
grid on;
xticks(0:0.2:1);
xticklabels({'0', '0.2\pi', '0.4\pi', '0.6\pi', '0.8\pi', '\pi'});


%% Part C(b): 2D Image Subband Decomposition with Perfect Reconstruction
% =============================================
%  Step 1: Load Image
%  =============================================

img = imread('owl.png'); 
img = im2double(img);  % Now in [0,1] range

%% =============================================
%  Step 2: Design a 2-Channel Perfect Reconstruction Filter Bank
%  =============================================
N = 9;  
band_edge = 0.45; % Cutoff frequency (normalized)
[h0, h1, g0, g1] = firpr2chfb(N, band_edge);

%% =============================================
%  Step 3: Image Sub-Band Decomposition
%  =============================================
% Apply filter bank row-wise
low_rows = conv2(img, h0(:), 'same');
high_rows = conv2(img, h1(:), 'same');

% Apply filter bank column-wise and downsample
LL = conv2(low_rows, h0(:)', 'same');  LL = LL(1:2:end, 1:2:end);
LH = conv2(low_rows, h1(:)', 'same');  LH = LH(1:2:end, 1:2:end);
HL = conv2(high_rows, h0(:)', 'same'); HL = HL(1:2:end, 1:2:end);
HH = conv2(high_rows, h1(:)', 'same'); HH = HH(1:2:end, 1:2:end);

% Display sub-band images
figure('Name', 'Subband Decomposition');
subplot(2,2,1); imshow(LL, []); title('LL (Low-Low)');
subplot(2,2,2); imshow(LH, []); title('LH (Low-High)');
subplot(2,2,3); imshow(HL, []); title('HL (High-Low)');
subplot(2,2,4); imshow(HH, []); title('HH (High-High)');

%% =============================================
%  Step 4: Exact Reconstruction from Sub-Bands
%  =============================================
% Upsample sub-bands
LL_up = zeros(512); LL_up(1:2:end, 1:2:end) = LL;
LH_up = zeros(512); LH_up(1:2:end, 1:2:end) = LH;
HL_up = zeros(512); HL_up(1:2:end, 1:2:end) = HL;
HH_up = zeros(512); HH_up(1:2:end, 1:2:end) = HH;

% Apply synthesis filters
recon_LL = conv2(conv2(LL_up, g0(:)', 'same'), g0(:), 'same');
recon_LH = conv2(conv2(LH_up, g1(:)', 'same'), g0(:), 'same');
recon_HL = conv2(conv2(HL_up, g0(:)', 'same'), g1(:), 'same');
recon_HH = conv2(conv2(HH_up, g1(:)', 'same'), g1(:), 'same');

% Combine reconstructed components
recon = recon_LL + recon_LH + recon_HL + recon_HH;

% Display reconstructed image
figure('Name', 'Reconstruction');
subplot(1,2,1); imshow(img); title('Original Image');
subplot(1,2,2); imshow(recon); title('Reconstructed Image');

% Compute SNR
error = img - recon;
SNR = 20*log10(norm(img(:))/norm(error(:)));
disp(['Reconstruction SNR: ', num2str(SNR), ' dB']);

%% =============================================
%  Step 5: Reconstruction with LH, HL, HH Set to Zero
%  =============================================
LL_only = recon_LL; % Only LL is used
figure('Name', 'LL Only Reconstruction');
imshow(LL_only); title('Reconstruction using only LL');

%% =============================================
%  Step 6: Quantization of LH, HL, HH Sub-Bands
%  =============================================
% Define quantization levels
quant_levels = {[-1, 1], [-1, 0, 1], [-1, -0.5, 0, 0.5, 1]}; % 2, 3, and 4 levels

for q = 1:length(quant_levels)
    levels = quant_levels{q};
    LH_q = quantize_subband(LH, levels);
    HL_q = quantize_subband(HL, levels);
    HH_q = quantize_subband(HH, levels);
    LH_up_q = zeros(512); LH_up_q(1:2:end, 1:2:end) = LH_q;
    HL_up_q = zeros(512); HL_up_q(1:2:end, 1:2:end) = HL_q;
    HH_up_q = zeros(512); HH_up_q(1:2:end, 1:2:end) = HH_q;
    recon_LH_q = conv2(conv2(LH_up_q, g1(:)', 'same'), g0(:), 'same');
    recon_HL_q = conv2(conv2(HL_up_q, g0(:)', 'same'), g1(:), 'same');
    recon_HH_q = conv2(conv2(HH_up_q, g1(:)', 'same'), g1(:), 'same');
    recon_q = recon_LL + recon_LH_q + recon_HL_q + recon_HH_q;

    figure('Name', ['Reconstruction with ', num2str(length(levels)), ' Quantization Levels']);
    imshow(recon_q);
    title(['Reconstruction with ', num2str(length(levels)), ' Quantization Levels']);
end

%% =============================================
%  Step 7: Further Decomposition of LL
%  =============================================
LL_low_rows = conv2(LL, h0(:), 'same');
LL_high_rows = conv2(LL, h1(:), 'same');
LL_LL = conv2(LL_low_rows, h0(:)', 'same');  LL_LL = LL_LL(1:2:end, 1:2:end);
LL_LH = conv2(LL_low_rows, h1(:)', 'same');  LL_LH = LL_LH(1:2:end, 1:2:end);
LL_HL = conv2(LL_high_rows, h0(:)', 'same'); LL_HL = LL_HL(1:2:end, 1:2:end);
LL_HH = conv2(LL_high_rows, h1(:)', 'same'); LL_HH = LL_HH(1:2:end, 1:2:end);

% Display second-level sub-bands
figure('Name', 'Second-Level Decomposition of LL');
subplot(2,2,1); imshow(LL_LL, []); title('LL-LL');
subplot(2,2,2); imshow(LL_LH, []); title('LL-LH');
subplot(2,2,3); imshow(LL_HL, []); title('LL-HL');
subplot(2,2,4); imshow(LL_HH, []); title('LL-HH');

disp('Processing Complete!');


% Extra code:
% %% Effect of Quantization on Second-Level Subbands
% clc; clear; close all;
% 
% % Load and prepare image
% img = im2double(imread('owl.png'));
% if size(img,3) == 3, img = rgb2gray(img); end
% img = imresize(img, [512 512]);
% 
% % Design filters
% N = 9;
% [h0, h1, g0, g1] = firpr2chfb(N, 0.45);
% 
% % First-level decomposition (just to get LL)
% LL = conv2(conv2(img, h0(:), 'same'), h0(:)');
% LL = LL(1:2:end, 1:2:end); % 256x256
% 
% % Second-level decomposition
% LL_LH = conv2(conv2(LL, h0(:), 'same'), h1(:)');
% LL_LH = LL_LH(1:2:end, 1:2:end); % 128x128
% 
% % Quantization levels to test
% quant_levels = [2 3 5]; 
% 
% % Create figure
% figure('Position', [100 100 1200 400]);
% 
% for i = 1:length(quant_levels)
%     % Quantize LL-LH subband
%     q_subband = mat2gray(LL_LH); % Normalize to [0,1]
%     q_subband = imquantize(q_subband, multithresh(q_subband, quant_levels(i)-1));
% 
%     % Plot
%     subplot(1, length(quant_levels), i);
%     imshow(q_subband, []);
%     title([num2str(quant_levels(i)) '-Level Quantization']);
% 
%     % Add PSNR calculation
%     orig = mat2gray(LL_LH);
%     quant = im2double(q_subband);
%     mse = mean((orig(:) - quant(:)).^2);
%     psnr = 10*log10(1/mse);
%     xlabel(['PSNR: ' num2str(psnr,2) ' dB']);
% end
% 
% % Add original for comparison
% figure;
% subplot(1,2,1); imshow(mat2gray(LL_LH)); title('Original LL-LH');
% subplot(1,2,2); histogram(LL_LH(:), 50); title('LL-LH Histogram');