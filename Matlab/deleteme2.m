%% Part (a) - Analytical DTFT expression
% x[n] = (0.7)^n for 0 <= n <= 7, 0 elsewhere

syms w n
a = 0.7;
Xw = symsum(a^n * exp(-1j*w*n), n, 0, 7);
Xw_mag = abs(Xw);
Xw_phase = angle(Xw);

% Plot magnitude and phase
f = linspace(-pi, pi, 1000);
Xw_eval = double(subs(Xw, w, f));

figure;
subplot(2,1,1);
plot(f, abs(Xw_eval));
title('DTFT Magnitude |X(w)|');
xlabel('\omega'); ylabel('|X(\omega)|'); grid on;

subplot(2,1,2);
plot(f, angle(Xw_eval));
title('DTFT Phase ∠X(w)');
xlabel('\omega'); ylabel('∠X(\omega)'); grid on;

%% Part (b) - 8-point DFT using fft
x = (0.7).^(0:7);  % original sequence of length 8
X8 = fft(x, 8);

figure;
subplot(2,1,1);
stem(0:7, abs(X8));
title('8-point DFT Magnitude');
xlabel('k'); ylabel('|X[k]|'); grid on;

subplot(2,1,2);
stem(0:7, angle(X8));
title('8-point DFT Phase');
xlabel('k'); ylabel('∠X[k]'); grid on;

%% Part (c) - 16-point DFT with zero-padding
x_pad16 = [x, zeros(1,8)];  % pad to length 16
X16 = fft(x_pad16);

figure;
subplot(2,1,1);
stem(0:15, abs(X16));
title('16-point DFT Magnitude (Zero-Padded)');
xlabel('k'); ylabel('|X[k]|'); grid on;

subplot(2,1,2);
stem(0:15, angle(X16));
title('16-point DFT Phase (Zero-Padded)');
xlabel('k'); ylabel('∠X[k]'); grid on;

% Comment: Zero-padding increases frequency resolution without changing content

%% Part (d) - 128-point DFT with zero-padding
x_pad128 = [x, zeros(1,120)];
X128 = fft(x_pad128);

figure;
subplot(2,1,1);
plot(0:127, abs(X128));
title('128-point DFT Magnitude');
xlabel('k'); ylabel('|X[k]|'); grid on;

subplot(2,1,2);
plot(0:127, angle(X128));
title('128-point DFT Phase');
xlabel('k'); ylabel('∠X[k]'); grid on;

%% Part (e) - Relationship between DTFT and DFT
% DFT samples are samples of the DTFT at uniform frequency intervals:
% ω_k = 2πk/N
% As N increases (e.g. from 8 to 128), DFT more closely approximates the DTFT


%% Part (a) - Load the signal and compute DFT with N = 25
s = load('tones.mat');
x = s.y1;

N1 = 25;
X1 = fft(x, N1);

figure;
plot(0:N1-1, abs(X1));
title('Magnitude of DFT (N = 25)');
xlabel('DFT Index k'); ylabel('|X[k]|'); grid on;

% Comment manually based on peaks: How many distinct frequency peaks are visible?

%% Part (b) - Zero-padding and comparison

% Try different zero-padded lengths
N2 = 50;
X2 = fft(x, N2);

figure;
plot(0:N2-1, abs(X2));
title('Magnitude of DFT (N = 50)');
xlabel('DFT Index k'); ylabel('|X[k]|'); grid on;

N3 = 100;
X3 = fft(x, N3);

figure;
plot(0:N3-1, abs(X3));
title('Magnitude of DFT (N = 100)');
xlabel('DFT Index k'); ylabel('|X[k]|'); grid on;

N4 = 256;
X4 = fft(x, N4);

figure;
plot(0:N4-1, abs(X4));
title('Magnitude of DFT (N = 256)');
xlabel('DFT Index k'); ylabel('|X[k]|'); grid on;

