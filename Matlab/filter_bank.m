function [y, y_no_delay] = filter_bank(s, g1, g2, g3, g4)
% FILTER_BANK Implements a 3-stage perfect reconstruction filter bank
%   [y, y_no_delay] = filter_bank(s, g1, g2, g3, g4) processes input signal s
%   through analysis and synthesis filter banks with gains g1-g4
%
%   Inputs:
%       s  - Input signal (row vector)
%       g1 - Gain for first synthesis stage low band
%       g2 - Gain for first synthesis stage high band
%       g3 - Gain for second synthesis stage high band
%       g4 - Gain for third synthesis stage high band
%
%   Outputs:
%       y          - Full reconstructed signal
%       y_no_delay - Reconstructed signal with delay compensation

    % Design filters (same as in main code)
    n = 21;
    [h_lp, h_hp, G_lp, G_hp] = firpr2chfb(n, 0.45);
    h = dfilt.delay(1000);  % Delay filter for phase compensation

    %% Analysis Stage (3 levels)
    % Stage 1
    s_low_1 = conv(s, h_lp, 'same');
    s_high_1 = conv(s, h_hp, 'same');
    sd_low_1 = s_low_1(1:2:end);  % Downsample
    sd_high_1 = s_high_1(1:2:end);

    % Stage 2
    s_low_2 = conv(sd_low_1, h_lp, 'same');
    s_high_2 = conv(sd_low_1, h_hp, 'same');
    sd_low_2 = s_low_2(1:2:end);  % Downsample
    sd_high_2 = s_high_2(1:2:end);

    % Stage 3
    s_low_3 = conv(sd_low_2, h_lp, 'same');
    s_high_3 = conv(sd_low_2, h_hp, 'same');
    sd_low_3 = s_low_3(1:2:end);  % Downsample
    sd_high_3 = s_high_3(1:2:end);

    %% Synthesis Stage (3 levels)
    % Stage 1
    su_low_3 = filter(h, sd_low_3) * g1;
    su_low_3(2,:) = 0;  % Upsample by 2
    su_low_3 = su_low_3(:).';
    
    su_high_3 = sd_high_3 * g2;
    su_high_3(2,:) = 0;  % Upsample by 2
    su_high_3 = su_high_3(:).';
    
    y_low_1 = conv(su_low_3, G_lp, 'same');
    y_high_1 = conv(su_high_3, G_hp, 'same');
    y_1 = y_high_1 + y_low_1;

    % Stage 2
    yu_1 = y_1;
    yu_1(2,:) = 0;  % Upsample by 2
    yu_1 = yu_1(:).';
    
    % Delay compensation for high band
    sd_high_2_delayed = [zeros(1,21), sd_high_2, zeros(1,21)];
    
    su_high_2 = filter(h, sd_high_2_delayed) * g3;
    su_high_2(2,:) = 0;  % Upsample by 2
    su_high_2 = su_high_2(:).';
    
    y_low_2 = conv(yu_1, G_lp, 'same');
    y_high_2 = conv(su_high_2, G_hp, 'same');
    y_2 = y_low_2 + y_high_2;

    % Stage 3
    yu_2 = y_2;
    yu_2(2,:) = 0;  % Upsample by 2
    yu_2 = yu_2(:).';
    
    % Delay compensation for high band
    sd_high_1_delayed = [zeros(1,63), sd_high_1, zeros(1,63)];
    
    su_high_1 = filter(h, sd_high_1_delayed) * g4;
    su_high_1(2,:) = 0;  % Upsample by 2
    su_high_1 = su_high_1(:).';
    
    y_low_3 = conv(yu_2, G_lp, 'same');
    y_high_3 = conv(su_high_1, G_hp, 'same');
    y = y_low_3 + y_high_3;

    % Delay-compensated output
    y_no_delay = y(148:147+length(s));
end