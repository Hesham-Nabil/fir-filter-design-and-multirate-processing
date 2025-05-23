 clc, clear , close all
 load nspeech2.mat
 s = nspeech2';
 clear nspeech2

  h = dfilt.delay(1000);
 g1=1;g2=2;g3=1;g4=1;

 n = 21;
 [h_lp,h_hp,G_lp,G_hp] = firpr2chfb(n,.45);
 %   Analysis Stage 1
 s_low_1  = conv(s, h_lp);
 s_high_1 = conv(s, h_hp);
 % Downsample
 sd_low_1  = s_low_1(1:2:end);
 sd_high_1 = s_high_1(1:2:end);
 %  Analysis Stage 2
 s_low_2  = conv(sd_low_1, h_lp);
 s_high_2 = conv(sd_low_1, h_hp);
 % Downsample
 sd_low_2  = s_low_2(1:2:end) ;
 sd_high_2 = s_high_2(1:2:end);
 % Analysis Stage 3
 s_low_3 = conv(sd_low_2, h_lp);
 s_high_3= conv(sd_low_2, h_hp);
 % Downsample
 sd_low_3  = s_low_3(1:2:end) ;
 sd_high_3 = s_high_3(1:2:end);
 %  Synthesis Stage 1
 % upsample
 su_low_3 =filter(h,sd_low_3)*g1;su_low_3(2,:) = 0; su_low_3 = su_low_3(:).';
 su_high_3=sd_high_3*g2;su_high_3(2,:) = 0; su_high_3 = su_high_3(:).';
 % interpolate
 y_low_1 = conv(su_low_3,G_lp);
 y_high_1 = conv(su_high_3,G_hp);
 y_1 = y_high_1 + y_low_1;
 %  Synthesis Stage 2
 % upsample
 yu_1 = y_1;yu_1(2,:) = 0; yu_1 = yu_1(:).';
 % delay the high signal frequency by 21 samples
 sd_high_2_delayed =[zeros(1,21), sd_high_2, zeros(1,21)];
 % upsample
 su_high_2=filter(h,sd_high_2_delayed)*g3;su_high_2(2,:) = 0; 
 su_high_2 = su_high_2(:).';
 % interpolate
 y_low_2 = conv(yu_1,G_lp);
 y_high_2 = conv(su_high_2,G_hp);
 y_2 = y_low_2+y_high_2;
 % Synthesis Stage 3
 % upsample
 yu_2 = y_2;yu_2(2,:) = 0; yu_2 = yu_2(:).';
 % delay the high signal frequency by 63 samples
 sd_high_1_delayed =[zeros(1,63), sd_high_1, zeros(1,63)];
 % upsample
 su_high_1=filter(h,sd_high_1_delayed)*g4;su_high_1(2,:) = 0; 
su_high_1 = su_high_1(:).';
 % interpolate
 y_low_3 = conv(yu_2,G_lp);
 y_high_3 = conv(su_high_1,G_hp);
 y = y_low_3+y_high_3;
 y_no_delay = y(148:147+length(s));


N = 1024; 
[sf, w] = DTFT(s, N);     
[yf, w] = DTFT(y_no_delay, N);  

% Plotting
figure(1)
subplot(3,2,1)
plot(s)
title('Input signal in time domain')
xlabel('n'); ylabel('signal')

subplot(3,2,2)
plot(y_no_delay)
title('Output signal in time domain')
xlabel('n'); ylabel('signal')


subplot(3,2,3)
plot(w, abs(sf))
title('Input signal amplitude Response')
xlabel('w (rad/s)'); ylabel('Amplitude Response')

subplot(3,2,4)
plot(w, abs(sf))
title('Input signal amplitude Response zoomed in')
xlabel('w (rad/s)'); ylabel('Amplitude Response')
xlim([0, pi])
xline(1/8*pi); xline(1/4*pi); xline(1/2*pi);


subplot(3,2,5)
plot(w, abs(yf))
title('Output signal amplitude Response')
xlabel('w (rad/s)'); ylabel('Amplitude Response')
xline(1/8*pi); xline(1/4*pi); xline(1/2*pi);


subplot(3,2,6)
plot(w, abs(yf))
title('Output signal amplitude Response Zoomed in')
xlabel('w (rad/s)'); ylabel('Amplitude Response')
xline(1/8*pi); xline(1/4*pi); xline(1/2*pi);
xlim([0, pi])

figure(2)
subplot(3,1,1)
plot(0:length(y_low_1)-1, y_low_1)
hold on
plot(0:length(y_high_1)-1, y_high_1)
title('Output of Synthesis level 3 (Delay=500 samples)')
xlabel('n'); ylabel('signal');
legend('Low freq', 'High freq')

subplot(3,1,2)
plot(0:length(y_low_2)-1, y_low_2)
hold on
plot(0:length(y_high_2)-1, y_high_2)
title('Output of Synthesis level 2')
xlabel('n'); ylabel('signal');
legend('Low freq', 'High freq')

subplot(3,1,3)
plot(0:length(y_low_3)-1, y_low_3)
hold on
plot(0:length(y_high_3)-1, y_high_3)
title('Output of Synthesis level 3')
xlabel('n'); ylabel('signal');
legend('Low freq', 'High freq')

figure(3)
n = 0:49;
plot(s)
hold on
plot(y_no_delay)
legend('Input signal', 'Output signal')
xlabel('n'); ylabel('signal')
title('Input and output signals no delay')
