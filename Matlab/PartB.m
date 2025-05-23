n=21
[h_lp,h_hp,G_lp,G_hp] = firpr2chfb(n,.45);
N = 1024;  

[h0, w] = DTFT(h_lp, N);
[h1, w] = DTFT(h_hp, N);
[g0, w] = DTFT(G_lp, N);
[g1, w] = DTFT(G_hp, N);

 subplot(2,2,1)
 plot(w,abs(h0));xlabel('w');title('H0 magnitude response')
 ylabel('Magnitude response')
 subplot(2,2,2)
 plot(w,abs(h1));xlabel('w');title('H1 magnitude response')
 ylabel('Magnitude response')
 subplot(2,2,3)
 plot(w,abs(g0));xlabel('w');title('G0 magnitude response')
 ylabel('Magnitude response')
 subplot(2,2,4)
 plot(w,abs(g1));xlabel('w');title('G1 magnitude response')
 ylabel('Magnitude response')


 fs= 48000;
 f_sig1 = 3000;  % Band 1
 f_sig2 = 8000;  % Band 2
 f_sig3 = 18000; % Band 3
 f_sig4 = 32000; % Band 4
 n = 0:100;
 sig1 = 0.7*sin(2* pi* f_sig1/fs .*n); % Band 1
 sig2 = 0.5*cos(2* pi* f_sig2/fs .*n); % Band 2
 sig3 = 0.3*cos(2* pi* f_sig3/fs .*n); % Band 3
 sig4 = 0.2*cos(2* pi* f_sig4/fs .*n); % Band 4
 s = sig1+sig2+sig3+sig4;
 [Sf, w]  = DTFT(s,N);
 subplot(2,1,1)
 stem(n,s); xlabel('n'); ylabel('signal');
 title('Artifitial signal of 4 frequencies')
 subplot(2,1,2)
 plot(w/pi,abs(Sf)); xlabel('w (rad/s)'); ylabel('Amplitude');
 title('Amplitude of the signal in the frequency domain')

 s_low_1  = conv(s, h_lp);
 s_high_1 = conv(s, h_hp);
 % Downsample
 sd_low_1  = s_low_1(1:2:end);
 sd_high_1 = s_high_1(1:2:end);
 s_low_2  = conv(sd_low_1, h_lp);
 s_high_2 = conv(sd_low_1, h_hp);
 % Downsample
 sd_low_2  = s_low_2(1:2:end) ;
 sd_high_2 = s_high_2(1:2:end);
 s_low_3 = conv(sd_low_2, h_lp);
 s_high_3= conv(sd_low_2, h_hp);
 % Downsample
 sd_low_3  = s_low_3(1:2:end) ;
 sd_high_3 = s_high_3(1:2:end);

 figure(1)
 subplot(3,1,1)
 plot(0:length(sd_low_1)-1, sd_low_1)
 hold on
 plot(0:length(sd_high_1)-1, sd_high_1)
 title('Output of Analysis level 1')
 xlabel('n'); ylabel('signal');legend('Low freq', 'High freq')
 subplot(3,1,2)
 plot(0:length(sd_low_2)-1, sd_low_2)
 hold on
 plot(0:length(sd_high_2)-1, sd_high_2)
 title('Output of Analysis level 2')
 xlabel('n'); ylabel('signal');legend('Low freq', 'High freq')
 subplot(3,1,3)
 plot(0:length(sd_low_3)-1, sd_low_3)
 hold on
 plot(0:length(sd_high_3)-1, sd_high_3)
 title('Output of Analysis level 3')
 xlabel('n'); ylabel('signal');legend('Low freq', 'High freq')

 su_low_3 =sd_low_3;su_low_3(2,:) = 0; su_low_3 = su_low_3(:).';
 su_high_3=sd_high_3;su_high_3(2,:) = 0; su_high_3 = su_high_3(:).';
 % interpolate
 y_low_1 = conv(su_low_3,G_lp);
 y_high_1 = conv(su_high_3,G_hp);
 y_1 = y_high_1 + y_low_1;

  yu_1 = y_1;yu_1(2,:) = 0; yu_1 = yu_1(:).';
 % delay the high signal frequency by 21 samples
 sd_high_2_delayed =[zeros(1,21), sd_high_2, zeros(1,21)];
 % upsample
 su_high_2=sd_high_2_delayed;su_high_2(2,:) = 0; su_high_2 = su_high_2(:).';

 % interpolate
 y_low_2 = conv(yu_1,G_lp);
 y_high_2 = conv(su_high_2,G_hp);
 y_2 = y_low_2+y_high_2;

  yu_2 = y_2;yu_2(2,:) = 0; yu_2 = yu_2(:).';
 % delay the high signal frequency by 63 samples
 sd_high_1_delayed =[zeros(1,63), sd_high_1, zeros(1,63)];
 % upsample
 su_high_1=sd_high_1_delayed;su_high_1(2,:) = 0; su_high_1 = su_high_1(:).';
 % interpolate
 y_low_3 = conv(yu_2,G_lp);
 y_high_3 = conv(su_high_1,G_hp);
 y = y_low_3+y_high_3;

  figure(2)
 subplot(3,1,1)
 plot(0:length(y_low_1)-1, y_low_1)
 hold on
 plot(0:length(y_high_1)-1, y_high_1)
 title('Output of Synthesis level 3')
 xlabel('n'); ylabel('signal');legend('Low freq', 'High freq')
 subplot(3,1,2)
 plot(0:length(y_low_2)-1, y_low_2)
 hold on
 plot(0:length(y_high_2)-1, y_high_2)
 title('Output of Synthesis level 2')
 xlabel('n'); ylabel('signal');legend('Low freq', 'High freq')
 subplot(3,1,3)
 plot(0:length(y_low_3)-1, y_low_3)
 hold on
 plot(0:length(y_high_3)-1, y_high_3)
 title('Output of Synthesis level 3')
 xlabel('n'); ylabel('signal');legend('Low freq', 'High freq')


  figure(3)
 n = 0:49;
 plot(s)
 hold on
 plot(y)
 legend('Input signal', 'Output signal')
 xlabel('n');ylabel('signal')
 title('Input and output signals')

 figure(4)
 subplot(2,1,1)
 [sf,w] = DTFT(s,N);
 plot(w,abs(sf));xlabel('w (rad/s)'); ylabel('Amplitude Response')
 title('Amplitude Response of the input signal')
 subplot(2,1,2)
 [yf,w] = DTFT(y,N);
 plot(w,abs(yf));xlabel('w (rad/s)'); ylabel('Amplitude Response')
 title('Amplitude Response of the output signal')

  figure(5)
 subplot(2,1,1)
 [sf,w] = DTFT(s_high_1,N);
 plot(w,abs(sf));xlabel('w (rad/s)'); ylabel('Amplitude Response')
 title('Amplitude Response of the high frequency signal in level 1 before downsampling')
 subplot(2,1,2)
 [yf,w] = DTFT(sd_high_1,N);
 plot(w,abs(yf));xlabel('w (rad/s)'); ylabel('Amplitude Response')
 title('Amplitude Response of the high frequency signal in level 1 after downsampling')



 

% Define gain factors for two equalizer settings

% Equalizer 1: Bass Boost (Enhance low frequencies, reduce highs)
G1_bass = 1.8;  
G2_bass = 1.2;  
G3_bass = 0.8;  
G4_bass = 0.5;  

% Equalizer 2: Treble Boost (Reduce lows, enhance highs)
G1_treble = 0.7;  
G2_treble = 0.9;  
G3_treble = 1.2;  
G4_treble = 1.8;  

% Apply gain factors for Bass Boost
sd_low_3_bass  = G1_bass * sd_low_3;
sd_high_3_bass = G2_bass * sd_high_3;
sd_high_2_bass = G3_bass * sd_high_2;
sd_high_1_bass = G4_bass * sd_high_1;

% Apply gain factors for Treble Boost
sd_low_3_treble  = G1_treble * sd_low_3;
sd_high_3_treble = G2_treble * sd_high_3;
sd_high_2_treble = G3_treble * sd_high_2;
sd_high_1_treble = G4_treble * sd_high_1;

% Function to match the lengths of two signals before addition
match_length = @(x, y) x(1:min(length(x), length(y))) + y(1:min(length(x), length(y)));

% Function for signal reconstruction
reconstruct = @(sd_low_3, sd_high_3, sd_high_2, sd_high_1) ...
    match_length( ...
        match_length( ...
            match_length(conv(sd_low_3, G_lp), conv(sd_high_3, G_hp)), ...
            conv(sd_high_2, G_hp)), ...
        conv(sd_high_1, G_hp) ...
    );

% Reconstruct equalized signals
y_eq_bass = reconstruct(sd_low_3_bass, sd_high_3_bass, sd_high_2_bass, sd_high_1_bass);
y_eq_treble = reconstruct(sd_low_3_treble, sd_high_3_treble, sd_high_2_treble, sd_high_1_treble);

% Ensure all signals have the same length before plotting
min_len = min([length(s), length(y_eq_bass), length(y_eq_treble)]);
s = s(1:min_len);
y_eq_bass = y_eq_bass(1:min_len);
y_eq_treble = y_eq_treble(1:min_len);

% Plot time-domain comparison
figure;
subplot(2,1,1)
plot(s, 'k', 'LineWidth', 1.2); hold on;
plot(y_eq_bass, 'b', 'LineWidth', 1.2);
plot(y_eq_treble, 'r', 'LineWidth', 1.2);
legend('Original', 'Bass Boost', 'Treble Boost');
xlabel('n'); ylabel('Signal Amplitude');
title('Time-Domain Comparison of Equalizers');
grid on;

% Frequency domain comparison
subplot(2,1,2)
[sf, w] = DTFT(s, N);
[yf_bass, ~] = DTFT(y_eq_bass, N);
[yf_treble, ~] = DTFT(y_eq_treble, N);
plot(w, abs(sf), 'k', 'LineWidth', 1.2); hold on;
plot(w, abs(yf_bass), 'b', 'LineWidth', 1.2);
plot(w, abs(yf_treble), 'r', 'LineWidth', 1.2);
legend('Original Spectrum', 'Bass Boosted Spectrum', 'Treble Boosted Spectrum');
xlabel('w (rad/s)'); ylabel('Magnitude');
title('Frequency-Domain Comparison of Equalizers');
grid on;





% Find the minimum length of the signals
min_len = min([length(s), length(y_eq_bass), length(y_eq_treble)]);

% Plot the first 100 ms of the signals
time_axis = (0:min_len-1) / fs; % Assuming fs is the sampling frequency (replace with the actual value)

% Get the first 100 ms of each signal (or the minimum length)
time_limit = round(0.1 * fs); % 100 ms in samples
if time_limit > min_len
    time_limit = min_len;
end

figure;
plot(time_axis(1:time_limit), s(1:time_limit), 'k', 'LineWidth', 1.2); hold on; % Original signal
plot(time_axis(1:time_limit), y_eq_bass(1:time_limit), 'b', 'LineWidth', 1.2); % Bass Boosted
plot(time_axis(1:time_limit), y_eq_treble(1:time_limit), 'r', 'LineWidth', 1.2); % Treble Boosted
legend('Original', 'Bass Boost', 'Treble Boost');
xlabel('Time [s]');
ylabel('Signal Amplitude');
title('Equalized Signals in Time Domain (First 100 ms)');
grid on;




%% additional plots
% Assuming n = 21 for the filter length, but we need to use the same length for the impulse response
n = 21;  % Number of samples for visualization
N = 1024;  % FFT length used for DTFT

% Compute the impulse responses of the decomposition filters (h_lp, h_hp)
impulse_h0 = ifft(h0, N); % Impulse response of h_lp (low-pass decomposition filter)
impulse_h1 = ifft(h1, N); % Impulse response of h_hp (high-pass decomposition filter)

% Compute the impulse responses of the synthesis filters (G_lp, G_hp)
impulse_g0 = ifft(g0, N); % Impulse response of G_lp (low-pass synthesis filter)
impulse_g1 = ifft(g1, N); % Impulse response of G_hp (high-pass synthesis filter)

% We limit the impulse responses to the first n samples for visualization
impulse_h0 = real(impulse_h0(1:n)); % Taking first n samples for visualization
impulse_h1 = real(impulse_h1(1:n));
impulse_g0 = real(impulse_g0(1:n));
impulse_g1 = real(impulse_g1(1:n));

% Ensure the time vector matches the length of the impulse response
time_vector = 0:n-1;  % Time vector for the first 'n' samples

% Ensure that time_vector and impulse responses have the same lengths
if length(impulse_h0) ~= length(time_vector)
    error('Length of time_vector and impulse_h0 do not match!');
end

% Plot the impulse responses
figure;

% Plot Impulse Response of H0 (Low-pass decomposition filter)
subplot(2,2,1)
stem(time_vector, impulse_h0, 'filled');
xlabel('Samples');
ylabel('Amplitude');
title('Impulse Response of H0 (Low-pass Decomposition)');

% Plot Impulse Response of H1 (High-pass decomposition filter)
subplot(2,2,2)
stem(time_vector, impulse_h1, 'filled');
xlabel('Samples');
ylabel('Amplitude');
title('Impulse Response of H1 (High-pass Decomposition)');

% Plot Impulse Response of G0 (Low-pass synthesis filter)
subplot(2,2,3)
stem(time_vector, impulse_g0, 'filled');
xlabel('Samples');
ylabel('Amplitude');
title('Impulse Response of G0 (Low-pass Synthesis)');

% Plot Impulse Response of G1 (High-pass synthesis filter)
subplot(2,2,4)
stem(time_vector, impulse_g1, 'filled');
xlabel('Samples');
ylabel('Amplitude');
title('Impulse Response of G1 (High-pass Synthesis)');
