function plot_filter_analysis(b, a, wp, ws, fs, noisy, filter_name, filter_order)
    % PLOT_FILTER_ANALYSIS Comprehensive filter analysis and plotting
    % Inputs:
    %   b, a - filter coefficients
    %   wp, ws - passband and stopband frequencies (rad/sample)
    %   fs - sampling frequency (Hz)
    %   noisy - input signal
    %   filter_name - string for title
    %   filter_order - filter order
    
    % Calculate multiplication operations
    mult_per_sample = 2 * filter_order + 1;
    fprintf('\n%s Filter:\n', filter_name);
    fprintf('Multiplication operations per sample (Direct Form II): %d\n', mult_per_sample);
    
    % Filter the signal
    filtered_sig = filter(b, a, noisy);
    
    % Create main analysis figure
    figure('Position', [100, 100, 900, 1000], 'Name', [filter_name ' Analysis']);
    
    % Magnitude response (dB)
    subplot(4,2,[1 2]);
    [H, w] = freqz(b, a, 1024);
    plot(w/pi, 20*log10(abs(H)), 'LineWidth', 1.5);
    title([filter_name ' (N=' num2str(filter_order) ') - Frequency Response']);
    xlabel('Normalized Frequency (\times\pi rad/sample)');
    xlim([0 1]);
    ylabel('Magnitude (dB)');
    grid on;
    hold on;

    % Mark passband and stopband edges
    passband_x = wp/pi;
    stopband_x = ws/pi;
    passband_y = 20*log10(abs(H(round(passband_x * 1024))));
    stopband_y = 20*log10(abs(H(round(stopband_x * 1024))));

    % Plot passband and stopband markers
    plot(passband_x, passband_y, 'ro', 'MarkerFaceColor', 'r');
    plot(stopband_x, stopband_y, 'ro', 'MarkerFaceColor', 'r');
    xline(passband_x, 'r--', 'Passband edge', 'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'middle', 'FontSize', 10);
    xline(stopband_x, 'r--', 'Stopband edge', 'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'middle', 'FontSize', 10);

    % Add horizontal text labels near markers
    text(passband_x - 0.04, passband_y - 10, 'Passband', 'Color', 'red', 'HorizontalAlignment', 'right');
    text(stopband_x + 0.06, stopband_y - 10, 'Stopband', 'Color', 'red', 'HorizontalAlignment', 'left');

    % Add legend for magnitude response
    legend('Magnitude Response', 'Passband Marker', 'Stopband Marker', 'Location', 'southwest');
    ylim([-100 5]);
    
    % Passband ripple (linear scale)
    subplot(4,2,3);
    plot(w/pi, abs(H), 'LineWidth', 1.5);
    xlim([0 passband_x * 1.1]);
    ylim([0.8 1.2]); % Y-axis from 0.8 to 1.2
    title('Passband Detail (Linear)');
    xlabel('Normalized Frequency');
    ylabel('Magnitude');
    grid on;
    hold on;
    plot(passband_x, abs(H(round(passband_x * 1024))), 'ro', 'MarkerFaceColor', 'r');
    text(passband_x, abs(H(round(passband_x * 1024))) + 0.02, 'Passband', 'HorizontalAlignment', 'right', 'Color', 'red');
    
    % Group delay plot
    subplot(4,2,4);
    [gd, w_gd] = grpdelay(b, a, 1024);
    plot(w_gd/pi, gd, 'LineWidth', 1.5);
    title('Group Delay');
    xlabel('Normalized Frequency');
    ylabel('Samples');
    grid on;
    xlim([0 1]);

    % Time domain comparison
    t = (0:length(noisy)-1)/fs;
    subplot(4,2,[5 6]);
    plot(t, noisy, 'b'); hold on;
    plot(t, filtered_sig, 'r', 'LineWidth', 1.2);
    title('Time Domain Comparison');
    xlabel('Time (s)');
    ylabel('Amplitude');
    xlim([0 0.5]); % First 0.5 seconds
    legend('Original', 'Filtered');
    grid on;
    
    % Frequency domain comparison
    Nfft = 4096;
    f = (0:Nfft/2-1)*fs/Nfft;
    P_orig = 20*log10(abs(fft(noisy,Nfft)));
    P_filt = 20*log10(abs(fft(filtered_sig,Nfft)));
    
    subplot(4,2,[7 8]);
    plot(f, P_orig(1:Nfft/2), 'b'); hold on;
    plot(f, P_filt(1:Nfft/2), 'r', 'LineWidth', 1.5);
    title('Frequency Domain Comparison');
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (dB)');
    legend('Original', 'Filtered');
    xlim([0 8000]);
    grid on;
    
end