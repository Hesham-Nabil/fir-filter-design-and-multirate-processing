function [X, w] = DTFT(x, N)
    % Compute the DTFT of a signal x at N frequency points
    % X: DTFT values
    % w: frequency points
    
    % Pad with zeros if necessary
    if length(x) < N
        x = [x(:); zeros(N-length(x), 1)];
    end
    
    % Compute DTFT using FFT
    X = fft(x, N);
    
    % Normalize
    X = X / length(x);
    
    % Frequency vector
    w = (0:N-1)/N * 2 * pi;
    w = w - 2*pi*(w>=pi); % Make it symmetric around 0
    
    % Shift to center 0 frequency
    X = fftshift(X);
    w = fftshift(w);
end