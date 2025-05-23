function plot_filter_response(H, w, wp, ws, filter_name)
    figure('Position', [100, 100, 900, 700]);
    
    % Magnitude Response (dB)
    subplot(3,1,1);
    plot(w/pi, 20*log10(abs(H)), 'LineWidth', 2);
    title([filter_name ' - Magnitude Response'], 'FontSize', 12);
    xlabel('Normalized Frequency (\times\pi rad/sample)', 'FontSize', 10);
    ylabel('Magnitude (dB)', 'FontSize', 10);
    grid on;
    hold on;
    line([wp/pi wp/pi], ylim, 'Color', 'r', 'LineStyle', '--');
    line([ws/pi ws/pi], ylim, 'Color', 'r', 'LineStyle', '--');
    xlim([0 1]);
    
    % Passband Ripple Detail (Linear)
    subplot(3,1,2);
    plot(w/pi, abs(H), 'LineWidth', 2);
    title('Passband Ripple Detail', 'FontSize', 12);
    xlabel('Normalized Frequency (\times\pi rad/sample)', 'FontSize', 10);
    ylabel('Magnitude', 'FontSize', 10);
    grid on;
    xlim([0 wp/pi*1.1]);
    
    % Group Delay
    subplot(3,1,3);
    [gd, w_gd] = grpdelay(H, 1, 512);
    plot(w_gd/pi, gd, 'LineWidth', 2);
    title('Group Delay', 'FontSize', 12);
    xlabel('Normalized Frequency (\times\pi rad/sample)', 'FontSize', 10);
    ylabel('Samples', 'FontSize', 10);
    grid on;
    xlim([0 1]);
end