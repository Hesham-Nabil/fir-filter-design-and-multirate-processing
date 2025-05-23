function h = LPFtrunc(N, wc)
    % Compute truncated impulse response of low-pass filter (centered)
    % N: filter size (must be odd)
    % wc: cutoff frequency
    
    M = (N-1)/2; % Middle point
    n = -M:M;    % Symmetric time indices around 0
    
    % Ideal impulse response (centered)
    hideal = (wc/pi) * sinc(wc/pi * n);
    
    % Rectangular window (simple truncation)
    h = hideal;
end