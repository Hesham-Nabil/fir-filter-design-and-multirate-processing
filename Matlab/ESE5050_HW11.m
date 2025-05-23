% Problem 1: LQR Design for Two-Trailer Robot

% Define system matrices
A = [0 0 0 0;
     1 -1 0 0;
     0 1 -1 0;
     1 0 0 0];

B = [1; 0; 0; 0];

% Choose diagonal Q (can adjust for state weighting)
Q = diag([10, 1, 1, 10]);

% Range of kappa values (R = scalar)
kappa_vals = logspace(-2, 2, 100);  % from 0.01 to 100
poles_all = zeros(length(kappa_vals), 4);

% Loop over each kappa value
for i = 1:length(kappa_vals)
    R = kappa_vals(i);
    K = lqr(A, B, Q, R);
    poles_all(i, :) = eig(A - B * K).';
end

% Plotting real vs. imaginary parts of poles (zoomed in)
figure;
hold on;
colors = lines(4);
for i = 1:4
    plot(real(poles_all(:, i)), imag(poles_all(:, i)), 'LineWidth', 2.5, 'Color', colors(i, :));
end
xlabel('Real Part', 'FontSize', 14);
ylabel('Imaginary Part', 'FontSize', 14);
title('Closed-Loop Pole Trajectories as \kappa Varies (Zoomed)', 'FontSize', 16);
xlim([-10 0]);  % Zoom into real axis from -10 to 0
ylim([-6 6]);
grid on;
legend('Pole 1', 'Pole 2', 'Pole 3', 'Pole 4', 'FontSize', 12);
set(gca, 'FontSize', 12);
line([0 0], ylim, 'Color', 'k', 'LineStyle', '-', 'LineWidth', 1);            % Imaginary axis (black)
line(xlim, [0 0], 'Color', [0.5 0.5 0.5], 'LineStyle', '--', 'LineWidth', 1); % Real axis (gray)

%% Problem 2
A = [0 0 0;
     0 0 1;
     0 -100 -0.2];
B = [1;
     0;
     0.1];
C = [1 0 -1];
D = 0;
[num, den] = ss2tf(A, B, C, D);
disp('Numerator:');
disp(num);
disp('Denominator:');
disp(den);
Gp = tf(num, den);
disp('Transfer Function Gp(s):');
Gp

%%

A = [0 0 0;
     0 0 1;
     0 -100 -0.2];
B = [1;
     0;
     0.1];

C1 = [B A*B A*A*B];
C2 = ctrb(A, B);
disp('Difference between C1 and C2:');
disp(C1 - C2);
disp('Rank of controllability matrix:');
rank(C1)
p = [-5, -0.1 + 10j, -0.1 - 10j];
K = place(A, B, p);
eig(A - B * K)

%%
% Define system matrices from Problem 2
A = [0 0 0;
     0 0 1;
     0 -100 -0.2];

B = [1;
     0;
     0.1];

C = [1 0 -1];
D = 0;

O = obsv(A, C);
rank_O = rank(O);
O_dual = ctrb(A', C')';
rank_O_dual = rank(O_dual);
disp(['Rank of observability matrix (obsv): ', num2str(rank_O)]);
disp(['Rank using duality: ', num2str(rank_O_dual)]);


p_obs = [-19, -20, -21];

L = place(A', C', p_obs)';
disp('Observer gain L =');
disp(L);

disp('Observer poles (eig(A - L*C)) =');
disp(eig(A - L * C));



% Original controller poles from assignment
p_K = [-5, -0.1 + 10j, -0.1 - 10j];

K = place(A, B, p_K);
disp('State-feedback gain K =');
disp(K);

disp('Controller poles (eig(A - B*K)) =');
disp(eig(A - B * K));


% Compensator dynamics
Ac = A - B*K - L*C;
Bc = L;
Cc = K;
Dc = 0;

[numc, denc] = ss2tf(Ac, Bc, Cc, Dc);
Gc = tf(numc, denc);

disp('Poles of Gc(s):');
disp(eig(Ac));


% Classical notch filter design
KGN = 5 * tf([1 0.2 100], [1 2 100]);

figure;
bode(KGN, Gc);
legend('Classical Notch', 'State-Space Compensator');
title('Bode Plot: Classical vs State-Space Compensator');
grid on;

% Recalculate K to increase damping (trial and error)
% Try poles with 10% damping (Zeta = 0.1) at 10 rad/s
Zeta = 0.1;
omega = 10;
p_revised = [-5, -Zeta*omega + 1i*omega*sqrt(1 - Zeta^2), -Zeta*omega - 1i*omega*sqrt(1 - Zeta^2)];

% New K
K_new = place(A, B, p_revised);
Ac_new = A - B*K_new - L*C;

% New compensator
[numc_new, denc_new] = ss2tf(Ac_new, Bc, K_new, 0);
Gc_new = tf(numc_new, denc_new);

% Display poles
disp('Poles of revised Gc(s):');
disp(eig(Ac_new));


figure;
bode(KGN, Gc, Gc_new);
legend('Classical Notch', 'Original State-Space', 'Revised State-Space');
title('Bode Plot Comparison: All Compensators');
grid on;




% Open-loop plant Gp(s)
Gp = tf([0.9 0.2 100], [1 0.2 100 0]);  % from earlier

% Loop gain for each controller
L_classical = KGN * Gp;
L_original = Gc * Gp;
L_revised = Gc_new * Gp;

% Loop bode plots with margins
figure;
margin(L_classical);
title('Loop Gain: Classical Notch');

figure;
margin(L_original);
title('Loop Gain: Original State-Space');

figure;
margin(L_revised);
title('Loop Gain: Revised State-Space');



% -------------------------------------------
% LOOP GAIN ANALYSIS USING CONTROLLER ONLY
% -------------------------------------------

% Open-loop plant
Gp = tf([0.9 0.2 100], [1 0.2 100 0]);

% Controller-only (no observer): Original K
Gc_open = ss(A - B*K, B, K, 0);
L_open = Gc_open * Gp;

% Controller-only: Revised K
Gc_open_new = ss(A - B*K_new, B, K_new, 0);
L_open_new = Gc_open_new * Gp;

% Plot margin for classical notch (as before)
figure;
margin(KGN * Gp);
title('Loop Gain: Classical Notch');

% Margin for original state-space controller without observer
figure;
margin(L_open);
title('Loop Gain: Original State-Space (No Observer)');

% Margin for revised state-space controller without observer
figure;
margin(L_open_new);
title('Loop Gain: Revised State-Space (10% Damping, No Observer)');
