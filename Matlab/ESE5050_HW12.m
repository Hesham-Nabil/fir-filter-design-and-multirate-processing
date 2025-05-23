%% Boeing 747 Yaw Damper System Analysis

% State-space model from assignment
A = [-0.05 -1.00 0.00 0.04;
      0.84 -0.15 -0.01 0.00;
     -3.00 0.41 -0.43 0.00;
      0.00 0.00 1.00 0.00];

B = [-0.01; 0.46; -0.11; 0.0];
C = [0 1 0 0];  % Output is yaw rate
D = 0;

% Aircraft plant
G = tf(ss(A, B, C, D));

% Washout filter
Hw = tf([2.7 0], [2.7 1]);

% Low-pass filter
Hlpf = tf([1], [0.27 1]);

% Open-loop system with filters
L = Hlpf * Hw * G;

% Now plot root locus for the real system
rlocus(L)
grid on;

rlocfind(L)


%%
% Aircraft plant (state-space to transfer function)
G = tf(ss(A, B, C, D));

% Washout filter
Hw = tf([2.7 0], [2.7 1]);

% Low-pass filter
Hlpf = tf([1], [0.27 1]);

% Open-loop transfer function (with filters)
L = Hlpf * Hw * G;

% Apply the proportional gain (K)
K = 1.2;  % Chosen gain
L_open_loop = K * L;

% Plot Bode plot with gain and phase margins
figure;
margin(L_open_loop);  % This also shows gain and phase margins
title(['Bode Plot of Open-Loop System with K = ' num2str(K)])
grid on;

% Optional: Display gain and phase margins
[GM, PM, Wcg, Wcp] = margin(L_open_loop);
fprintf('Gain Margin: %.2f dB at ω = %.2f rad/s\n', 20*log10(GM), Wcg);
fprintf('Phase Margin: %.2f° at ω = %.2f rad/s\n', PM, Wcp);

%%

% Define simulation time
T = 30;               % total time in seconds
N = 3001;             % number of samples

% Construct time vector from 0 to 30s, with 3001 points
t = linspace(0, T, N);

% input = out.simout;
% out_k_12 = out.simout1;

figure;
plot(t, out_k_12, 'b', 'LineWidth', 1.5); hold on;
plot(t, input, 'g', 'LineWidth', 1.5);

xlabel('Time (s)');
ylabel('Yaw Rate (rad/s)');
title('Yaw Rate Response at Different Gains');
legend('Response at K = 1.2', 'Input');
ylim([-2 2]);
grid on;

%% problem 2

A = [-0.05 -1.00 0.00 0.04;
      0.84 -0.15 -0.01 0.00;
     -3.00 0.41 -0.43 0.00;
      0.00 0.00 1.00 0.00];

eig_A = eig(A)


zeta = 0.7;
omega_n = sqrt((-0.0291)^2 + (0.9480)^2);

% New Dutch roll poles
real_part = -zeta * omega_n;
imag_part = omega_n * sqrt(1 - zeta^2);
new_dutch_roll = [real_part + 1i*imag_part, real_part - 1i*imag_part]


A = [-0.05 -1.00 0.00 0.04;
      0.84 -0.15 -0.01 0.00;
     -3.00  0.41 -0.43 0.00;
      0.00  0.00  1.00 0.00];

B = [-0.01; 0.46; -0.11; 0.0];

% Desired poles
p1 = -0.6639 + 0.6773i;
p2 = -0.6639 - 0.6773i;
p3 = -0.5634;
p4 = -0.0083;
desired_poles = [p1, p2, p3, p4];

% Compute state feedback gain K
K = place(A, B, desired_poles)


C = [0 1 0 0];  % Output is yaw rate
D = 0;

% Open-loop steady-state gain
G_ol = C * (-A)^(-1) * B;

% Closed-loop steady-state gain
A_cl = A - B * K;
G_cl = C * (-A_cl)^(-1) * B;

% Feedforward gain
K_delta = G_ol / G_cl

%%
% Given matrices
A = [-0.05 -1.00 0.00 0.04;
      0.84 -0.15 -0.01 0.00;
     -3.00  0.41 -0.43 0.00;
      0.00  0.00  1.00 0.00];

B = [-0.01; 0.46; -0.11; 0.0];
C = [0 1 0 0];
D = 0;

% Feedback and feedforward gains
K = [-0.2905, 2.7445, -0.0377, -0.1074];
K_delta = 0.9958;

% Closed-loop system
A_cl = A - B * K;
B_cl = B * K_delta;

% Define state-space system
sys_cl = ss(A_cl, B_cl, C, D);

% Simulate response to step or doublet input
t = linspace(0, 30, 3001);
u = zeros(size(t));           % zero until t = 1
u(t >= 1 & t < 4) = 1;        % +1 from 1s to 4s
u(t >= 4 & t < 7) = -1;       % -1 from 4s to 7s
u(t >= 7) = 0;                % return to zero

% Simulate response
[y, ~, ~] = lsim(sys_cl, u, t);

% Plot response
figure;
plot(t, y, 'b', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Yaw Rate (rad/s)');
title('Yaw Rate Response – State-Space Controller');
grid on;


figure;

% Plot classical response (K = 1.2)
plot(t, out_k_12, 'b', 'LineWidth', 1.5); hold on;

% Plot new state-space response
plot(t, y, 'r--', 'LineWidth', 1.5);  % Dashed red line

% Plot input (optional)
plot(t, input, 'g', 'LineWidth', 1.2);

xlabel('Time (s)');
ylabel('Yaw Rate (rad/s)');
title('Yaw Rate Comparison: Classical vs. State-Space Controller');
legend('Classical Design (K = 1.2)', 'State-Space Design', 'Input Signal');
ylim([-2 2]);
grid on;
