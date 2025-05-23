% Robotic Car with Trailers - Full State Feedback Control

clc; clear; close all;

%% Parameters and System Definition
a = 1.0;

A = [ 0   0   0   0;
      a  -a   0   0;
      0   a  -a   0;
      a   0   0   0 ];

B = [a; 0; 0; 0];
C = [0 -1 -1 1];
D = 0;

fast_poles = [-5, -2, -3 + 3j, -3 - 3j];
K_fast = place(A, B, fast_poles)
Kf_fast = 1 / (C * inv(A - B * K_fast) * B)


%% Confirm Controllability (Matrix + PBH)
ctrb_matrix = ctrb(A, B);
assert(rank(ctrb_matrix) == size(A, 1), 'System is not controllable using the controllability matrix.');

eigA = eig(A);
for i = 1:length(eigA)
    lambda = eigA(i);
    M = [lambda * eye(size(A)) - A, B];
    assert(rank(M) == size(A, 1), 'PBH test failed at eigenvalue %.2f', lambda);
end

%% Time and Input
Tmax = 30;
t = linspace(0, Tmax, 1000)';
y_d = ones(size(t));

%% Pole Set 1 (Slower response)
poles1 = [-0.2 + 0.3i, -0.2 - 0.3i, -0.9, -1.1];
K1 = place(A, B, poles1);
KF1 = 1 / (C * ((A - B*K1) \ B));

sys1 = ss(A - B*K1, B*KF1, C, D);
[y1, ~, x1] = lsim(sys1, y_d, t);
u1 = KF1 * y_d - x1 * K1';

%% Pole Set 2 (Faster response)
poles2 = [-5, -2, -3 + 3i, -3 - 3i];
K2 = place(A, B, poles2);
KF2 = 1 / (C * ((A - B*K2) \ B));

sys2 = ss(A - B*K2, B*KF2, C, D);
[y2, ~, x2] = lsim(sys2, y_d, t);
u2 = KF2 * y_d - x2 * K2';

%% Classical Controller (from last week)
Kc = 4;
Tc = 4;
G = tf(1, [1 2 1 0 0]);
Cc = Kc * tf([Tc 1], [1]);
G_open = Cc * G;
G_cl = feedback(G_open, 1);
[y_classical, t_classical] = step(G_cl, t);

%% Plot: Output Comparison
figure;
plot(t, y1, 'b', t, y2, 'r', t_classical, y_classical, 'k--', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Output y(t)');
title('Comparison of Output Responses');
legend('Pole Set 1', 'Pole Set 2', 'Classical Controller');
grid on;

%% Plot: Control Comparison
figure;
plot(t, u1, 'b', t, u2, 'r', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Control Input u(t)');
title('Comparison of Control Inputs');
legend('Pole Set 1', 'Pole Set 2');
grid on;

%% Root Locus and Bode Plot for Classical Controller
figure;
rlocus(G_open);

figure;
margin(G_open);   
[GM, PM, wcg, wcp] = margin(L);
grid on;



%% output hw8 vs output hw9
% slow part

Kf = 0.1287;
%output_signal_slow = out.output;
%output_classical = out.classical;
%input_signal_slow = out.input;
%input_classical = out.classicalinput;
T = 50;
dt = 0.01;
t = 0:dt:T;

figure
% --- Left: Output comparison ---
subplot(1, 2, 1)
plot(t, output_classical, 'b', 'LineWidth', 1.6)
hold on
plot(t, output_signal_slow, 'r', 'LineWidth', 1.6)
grid on
title('System Output y(t)')
xlabel('Time (s)')
ylabel('Output y')
legend('HW8: Classical', 'HW9: State-Feedback')
xlim([0, 40])

% --- Right: Input comparison ---
subplot(1, 2, 2)
plot(t, input_classical, 'b', 'LineWidth', 1.6)
hold on
plot(t, input_signal_slow, 'r', 'LineWidth', 1.6)
grid on
title('Control Input u(t)')
xlabel('Time (s)')
ylabel('Input u')
legend('HW8: Classical', 'HW9: State-Feedback')
ylim([-1.2, 1.2])
xlim([0, 10])


%% output hw8 vs output hw9 after optimizing:
% Output: HW8 vs HW9 (Fast Poles)
%output_signal_fast = out.output;
%input_signal_fast = out.input;

Kf = 180;
T = 50;
dt = 0.01;
t = 0:dt:T;

figure

% --- Left: Output comparison (fast poles) ---
subplot(1, 2, 1)
plot(t, output_classical, 'b', 'LineWidth', 1.6)
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


