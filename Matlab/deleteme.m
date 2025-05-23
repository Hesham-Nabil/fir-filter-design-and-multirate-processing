% Define system parameter
a = 1;

% Transfer function of plant: Gp(s) = a^4 / (s^2(s + a)^2)
num = a^4;
den = conv([1 0 0], conv([1 a], [1 a]));  % s^2(s+a)^2

Gp = tf(num, den);

% Define PD controller: C(s) = (Tc*s + 1)
Tc = 4;
s = tf('s');
C = (Tc*s + 1);

% Open-loop system
G_open = C * Gp;

% Plot root locus
figure
rlocus(G_open)
grid on
title('Root Locus of PD-Controlled System')

% Choose a gain K based on the plot
K = 0.09;  % example good gain

% Closed-loop system with unity feedback
G_cl = feedback(K * G_open, 1);

% Step response
figure
step(G_cl)
title('Step Response of Closed-Loop System (PD Control)')
ylabel('Output y(t)')
xlabel('Time (seconds)')
grid on

[GM, PM, Wcg, Wcp] = margin(G_open);

% Convert GM to dB
GM_dB = 20*log10(GM);

% Create title string with margin info
title_str = sprintf('Bode Plot (GM = %.2f dB @ %.2f rad/s, PM = %.2f° @ %.2f rad/s)', ...
    GM_dB, Wcg, PM, Wcp);

% Plot Bode with margin
figure
margin(G_open)
title(title_str)
grid on



%%%%%%%%%%%
T = 50;
dt = 0.01;
t = 0:dt:T;
u = ones(size(t));  % unit step input
% Closed-loop system
G_cl = feedback(K * C * Gp, 1);

% Generate step input
u = ones(size(t));

% Simulate system response using lsim on same time vector
y_pd = lsim(G_cl, u, t);

% Now plot it alongside your other signal
figure

% --- Left: Output y(t) comparison ---
subplot(1, 2, 1)
plot(t, y_pd, 'b', 'LineWidth', 1.5)
hold on
plot(t, output_signal_slow, 'r', 'LineWidth', 1.5)
grid on
title('System Output y(t)')
xlabel('Time (seconds)')
ylabel('Output y(t)')
legend('HW8: Classical PD','HW9: State-Feedback')
xlim([0, 50])
ylim([0, 2])

% --- Right: Control input u(t) comparison ---
subplot(1, 2, 2)
plot(t, input_classical, 'b', 'LineWidth', 1.6)
hold on
plot(t, input_signal_slow, 'r', 'LineWidth', 1.6)
grid on
title('Control Input u(t)')
xlabel('Time (s)')
ylabel('Control u(t)')
legend('HW8: Classical', 'HW9: State-Feedback')
ylim([-1.2, 1.2])
xlim([0, 10])




figure

% --- Left: Output comparison (fast poles) ---
subplot(1, 2, 1)
plot(t, y_pd, 'b', 'LineWidth', 1.5)
hold on
plot(t, output_signal_fast, 'r', 'LineWidth', 1.6)
grid on
title('System Output y(t)')
xlabel('Time (s)')
ylabel('Output y(t)')
legend('HW8: Classical', 'HW9: Fast State-Feedback')
xlim([0, 40])

% --- Right: Input comparison (fast poles) ---
subplot(1, 2, 2)
plot(t, input_classical, 'b', 'LineWidth', 1.6)
hold on
plot(t, input_signal_fast, 'r', 'LineWidth', 1.6)
grid on
title('Control Input u(t)')
xlabel('Time (s)')
ylabel('Control Input u(t)')
legend('HW8: Classical', 'HW9: Fast State-Feedback')
ylim([-3, 3])
xlim([0, 10])
