%% Problem 1: LQR and State Estimation
clc;
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');

k1 = 1; k2 = 1; k3 = 1; T  = 0.1;

A = [ 1     T     ;
     -k2*T  1-k1*T; ];

B = [ 0;  k3*T; ];

C = [ 1  0 ];

Q = diag([ 4  4 ]); R = 1;

[K, ~, e] = dlqr(A, B, Q, R, []);
disp(e);

p_re = 0.5; p_im = 0.03;
p  = [ p_re + p_im*1i  p_re - p_im*1i ];
Kf = place(A', C', p)';

%% Problem 1b): Simulation

N = 50;
x     = zeros(2, N+1);
x_hat = zeros(2, N+1);
y     = zeros(1, N);
u     = zeros(1, N);

x(:, 1)     = [ 5; 1; ];
x_hat(:, 1) = [ 6; 0; ];

for t = 1:N
    u(t) = -K * x_hat(:, t);
    y(t) =  C * x(:, t);
    x(:, t+1)     = A * x(:, t) + B * u(t);
    x_hat(:, t+1) = A * x_hat(:, t) + B * u(t) + Kf * (y(t) - C * x_hat(:, t));
end

time = 0:N;

figure;
subplot(2,1,1);
plot(time, x(1,:), 'b-', 'LineWidth', 1.5); hold on;
plot(time, x_hat(1,:), 'r--', 'LineWidth', 1.5);
title('State $x_1$ vs Estimate $\hat{x}_1$');
legend('True $x_1$', 'Estimated $\hat{x}_1$');
xlabel('$t$');
ylabel('$x_1$')
grid on;

subplot(2,1,2);
plot(time, x(2,:), 'b-', 'LineWidth', 1.5); hold on;
plot(time, x_hat(2,:), 'r--', 'LineWidth', 1.5);
title('State $x_2$ vs Estimate $\hat{x}_2$');
legend('True $x_2$', 'Estimated $\hat{x}_2$');
xlabel('$t$');
ylabel('$x_1$')
grid on;

%% Problemm 1c): LQG Eigenvalues

Phi = [ A - B*K      B*K     ; 
        zeros(2, 2)  A - Kf*C; ];

lqg_eigs = eig(Phi);
disp(lqg_eigs);

%% Problem 2: MPC and State Estimation
clc;

N_mpc = 30;    % prediction horizon
T_sim = 50;

u_min = -4; u_max = 4;
[nx, nu] = size(B);

x0 = x(:, 1); x_hat0 = x_hat(:, 1);

[G, Aeq, ~, lb, ub] = build_qp_matrices(A, B, x_hat0, Q, R, N_mpc, u_min, u_max);

[x, x_hat, u] = run_mpc_observer_simulation(A, B, C, Kf, N_mpc, Aeq, G, T_sim, x0, x_hat0, nx, nu, lb, ub);

time_x = 0:T_sim; time_u = 0:T_sim-1;

figure('Name', 'MPC with State Estimation', 'Position', [100, 100, 600, 800]);

subplot(3,1,1);
plot(time_x, x(1,:), 'b-o', 'LineWidth', 1.5); hold on;
plot(time_x, x_hat(1,:), 'r--', 'LineWidth', 1.5);
title('System State: $x_1$ vs Estimate $\hat{x}_1$');
ylabel('$x_1$');
legend('True $x_1$', 'Estimated $\hat{x}_1$');
grid on;

subplot(3,1,2);
plot(time_x, x(2,:), 'b-o', 'LineWidth', 1.5); hold on;
plot(time_x, x_hat(2,:), 'r--', 'LineWidth', 1.5);
title('System State: $x_2$ vs Estimate $\hat{x}_2$');
ylabel('$x_2$');
legend('True $x_2$', 'Estimated $\hat{x}_2$');
grid on;

subplot(3,1,3);
stairs(time_u, u, 'k-o', 'LineWidth', 1.5); hold on;
yline(u_max, 'm--');
yline(u_min, 'm--');
title('Control Input: $u$ ($-4 \leq u_t \leq 4$)');
ylabel('$u$');
xlabel('$t$');
ylim([-5 5]); 
grid on;

%% Problem 3: Infinite Horizon MPC
[~, P, ~]  = dlqr(A, B, Q, R);

N_mpc = 10;

[G, Aeq, ~, lb, ub] = build_qp_matrices_terminal(A, B, x0, Q, R, P, N_mpc, u_min, u_max);

[x, u] = run_mpc_simulation(A, B, N_mpc, Aeq, G, T_sim, x0, nx, nu, lb, ub);

figure('Name', 'Infinite-Horizon MPC', 'Position', [100, 100, 600, 600]);

subplot(2,1,1);
plot(time_x, x(1,:), 'b-o', 'LineWidth', 1.5); hold on;
plot(time_x, x(2,:), 'r-o', 'LineWidth', 1.5);
title('System States (Terminal Cost $P$, $N=10$)');
legend('$x_1$', '$x_2$');
xlabel('$t$');
ylabel('$x$');
grid on;

subplot(2,1,2);
stairs(time_u, u, 'k-o', 'LineWidth', 1.5); hold on;
yline(u_max, 'm--');
yline(u_min, 'm--');
title('Control Input: $u$');
ylim([-5 5]); 
xlabel('$t$');
ylabel('$u$')
grid on;





























