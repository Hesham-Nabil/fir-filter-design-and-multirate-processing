%% Problem 2: Equivalent Compensator Comparison
clear all;
% Define system matrices (same as HW09/10)
A = [0 0 0 0;
     1 -1 0 0;
     0 1 -1 0;
     1 0 0 0];

B = [1; 0; 0; 0];
C = [0 -1 -1 1];  % Output: distance of second trailer from wall


pc = [-1.1 -0.9 -0.5+0.5j -0.5-0.5j];
K = place(A, B, pc);

s = tf('s');

Gp = 1 / (s^2 * (s + 1)^2);  % Simple second-order system as an example

% Observer gain L from HW10
L = [15.8004; 14.8500; 10.9500; 31.8000];  % 4x1

% Equivalent compensator state-space matrices
Ac = A - B*K - L*C;
Bc = L;
Cc = K;
Dc = 0;

% Transfer function of equivalent compensator
[numc, denc] = ss2tf(Ac, Bc, Cc, Dc);
Gc = tf(numc, denc)

% Classical PD controller from HW8

%
% Design parameters
K2 = 0.09;        % Proportional gain
TC = 4;     % Time constant

s = tf('s');
Gc2 = K2 * (TC * s + 1)

Gp2 = 1 / (s^2 * (s + 1)^2);  % Simple second-order system as an example

C_pd = Gc2 * Gp2


% Bode plot of equivalent compensator and classical PD controller
figure;
bode(Gc, C_pd);
grid on;
legend('Equivalent Compensator', 'Classical PD Controller');
title('Bode Plot: Equivalent Compensator vs Classical PD Controller');
