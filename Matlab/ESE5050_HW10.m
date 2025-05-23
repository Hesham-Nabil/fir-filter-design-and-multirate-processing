A = [0 0 0 0;
     1 -1 0 0;
     0 1 -1 0;
     1 0 0 0];

B = [1; 0; 0; 0];

C = [0 -1 -1 1];

pc = [-1.1 -0.9 -0.5+0.5j -0.5-0.5j];

K = place(A, B, pc);

KF = 1 / (C * inv(A - B*K) * B);

O = obsv(A, C);

% Check rank
rank_O = rank(O);
disp(['Observability matrix rank: ', num2str(rank_O)]);


pe = [-2.2 -2.1 -1.9 -1.8];
L = place(A', C', pe)';

%% HW8 response
% Time vector
T = 50;
dt = 0.01;
t = (0:dt:T)';
n = length(t);

% PD Controller parameters (from your HW8)
K2 = 0.09;
Tc = 4;

% Define plant: G(s) = 1 / (s^2 * (s + 1)^2)
num_plant = 1;
den_plant = conv([1 0 0], conv([1 1], [1 1]));
G = tf(num_plant, den_plant);

% PD control is u(t) = K * (Tc * de/dt + e)

% Closed-loop system (using step for output)
C_pd = tf([K2*Tc, K2], [1]); % not used for input anymore
closed_loop = feedback(C_pd * G, 1);
[output_classical, ~] = step(closed_loop, t);

% Assign output
output_classical = output_classical(:);  % ensure column

% Compute error: e(t) = y_d - y(t)
e = 1 - output_classical;

% Compute derivative of error using finite differences
de_dt = [0; diff(e) / dt];

% Compute control input u(t) = K * (Tc * de/dt + e)
u = K2 * (Tc * de_dt + e);
input_classical = u(:);  % ensure column

% Plot output
figure;
subplot(2,1,1)
plot(t, output_classical, 'LineWidth', 1.5)
grid on
xlabel('Time (s)')
ylabel('Output y(t)')
title('HW8: Classical PD Controller Output')

% Plot control input
subplot(2,1,2)
plot(t, input_classical, 'LineWidth', 1.5)
grid on
xlabel('Time (s)')
ylabel('Control Input u(t)')
title('HW8: Classical PD Controller Input')







%% output hw8 vs output hw10

%output_signal_HW10 = out.output;
%input_signal_HW10 = out.input;

figure
% --- Left: Output comparison ---
subplot(1, 2, 1)
plot(t, output_classical, 'b', 'LineWidth', 1.6)
hold on
plot(t, output_signal_HW10, 'r', 'LineWidth', 1.6)
grid on
title('System Output y(t)')
xlabel('Time (s)')
ylabel('Output y')
legend('HW8: Classical', 'HW10: State-Feedback With estimator')
xlim([0, 50])

% --- Right: Input comparison ---
subplot(1, 2, 2)
plot(t, input_classical, 'b', 'LineWidth', 1.6)
hold on
plot(t, input_signal_HW10, 'r', 'LineWidth', 1.6)
grid on
title('Control Input u(t)')
xlabel('Time (s)')
ylabel('Input u')
legend('HW8: Classical', 'HW10: State-Feedback With estimator')
ylim([-1.2, 1.2])
xlim([0, 50])


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




%% Problem 2: Equivalent Compensator Comparison

% Define system matrices (same as HW09/10)
A = [0 0 0 0;
     1 -1 0 0;
     0 1 -1 0;
     1 0 0 0];

B = [1; 0; 0; 0];
C = [0 -1 -1 1];  % Output: distance of second trailer from wall

% Classical PD controller from HW8

%
% Design parameters
K2 = 0.09;        % Proportional gain
TC = 4;     % Time constant

s = tf('s');
Gc = K2 * (TC * s + 1);

Gp = 1 / (s^2 * (s + 1)^2);  % Simple second-order system as an example

C_pd = Gc * Gp;

% Feedback gain K from HW10 (pole placement)
K = [-1 0.005 -0.005 -0.495];  % 1x4

% Observer gain L from HW10
L = [15.8004; 14.8500; 10.9500; 31.8000];  % 4x1

% Equivalent compensator state-space matrices
Ac = A - B*K - L*C;
Bc = L;
Cc = K;
Dc = 0;

% Transfer function of equivalent compensator
[numc, denc] = ss2tf(Ac, Bc, Cc, Dc);
Gc = tf(numc, denc);

% Bode plot comparison
figure;
bode(C_pd, 'b', Gc, 'r')
grid on
legend('Classical PD: K(T_c s + 1)', 'Observer-Based: G_c(s)')
title('Bode Plot Comparison of Compensators')

% Compute poles of the compensator (internal dynamics)
eig_compensator = eig(Ac);
disp('Poles of A - BK - LC:')
disp(eig_compensator)
