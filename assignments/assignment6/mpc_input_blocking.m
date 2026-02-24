%% Problem 3: MPC and Input Blocking
close all; clc;

set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');

%% System Model
A = [ 0     0     0   ;
      0     0     1   ;
      0.1  -0.79  1.78;];

B = [ 1;  0;  0.1;];

C = [ 0  0  1];

x0 = [ 0;  0;  1;];

N = 30; r = 1;

Qt = 1; Rt = r;

nx = size(A,1); % nx: number of states
nu = size(B,2); % nu: number of controls
ny = size(C,1); % ny: number of controls

%% Problem 3a): Open-loop Optimization Problem
[G_a, Aeq_a, beq_a, lb_a, ub_a] = build_qp_matrices(A, B, C, x0, Qt, Rt, N);

z0_a = zeros(N*(nx+nu), 1); 

opt = optimset('Display','iter-detailed','Algorithm', 'active-set');
[z_a, ~, ~, output_a] = quadprog(G_a, [], [], [], Aeq_a, beq_a, lb_a, ub_a, z0_a, opt);
fprintf('Task (a) iterations: %d\n', output_a.iterations);

x_opt_a = reshape(z_a(1:N*nx), nx, N);
x_full_a = [x0, x_opt_a]; 
y_t_a = x_full_a(3, :);
    
u_opt_a = reshape(z_a(N*nx+1:end), nu, N);

plot_mpc_results(y_t_a, u_opt_a, N, 'Problem 3a: Open-Loop MPC');

%% Problem 3b): Input Blocking
nb = 6;
blocks_b = ones(1, nb) * 5;
[G_b, Aeq_b, beq_b, lb_b, ub_b, T_b] = build_qp_blocked(A, B, C, x0, Qt, Rt, blocks_b);

z0_b = zeros(N*nx + nb*nu, 1);

[z_b, ~, ~, output_b] = quadprog(G_b, [], [], [], Aeq_b, beq_b, lb_b, ub_b, z0_b, opt);
fprintf('Task (b) iterations: %d\n', output_b.iterations);

x_opt_b = reshape(z_b(1:N*nx), nx, N);
x_full_b = [x0, x_opt_b];
y_t_b = x_full_b(3, :);

u_blocked_b = T_b * z_b(N*nx+1:end);

plot_mpc_results(y_t_b, u_blocked_b, N, 'Problem 3b: Input Blocking')

%% Problem 3c): Increasing Input Block Lengths
blocks_c = [1, 1, 2, 4, 8, 14];
[G_c, Aeq_c, beq_c, lb_c, ub_c, T_c] = build_qp_blocked(A, B, C, x0, Qt, Rt, blocks_c);

z0_c = zeros(N*nx + nb*nu, 1);

[z_c, ~, ~, output_c] = quadprog(G_c, [], [], [], Aeq_c, beq_c, lb_c, ub_c, z0_c, opt);
fprintf('Task (c) iterations: %d\n', output_c.iterations);

x_opt_c = reshape(z_c(1:N*nx), nx, N);
x_full_c = [x0, x_opt_c];
y_t_c = x_full_c(3, :);

u_blocked_c = T_c * z_c(N*nx+1:end);

plot_mpc_results(y_t_c, u_blocked_c, N, 'Problem 3c: Input Blocking (Varying Lengths)')

%% Problem 3d): MPC Problem
[G_d, Aeq_d, beq_d, lb_d, ub_d] = build_qp_matrices(A, B, C, x0, Qt, Rt, N);

T_sim = 30;
[x_t_d, u_t_d] = run_mpc_simulation(A, B, N, Aeq_d, G_d, T_sim, x0, nx, nu, lb_d, ub_d);
y_t_d = x_t_d(3, :);

plot_mpc_results(y_t_d, u_t_d, N, 'Problem 3d: Closed-loop MPC Simulation')

%% Problem 3e): MPC With Input Blocking - Comparison
[x1_t_e, u1_t_e] = run_mpc_simulation(A, B, N, Aeq_b, G_b, T_sim, x0, nx, nu, lb_b, ub_b);
y1_t_e = x1_t_e(3, :);

plot_mpc_results(y1_t_e, u1_t_e, N, 'Problem 3e.1: MPC With Constant Input Blocks')

[x2_t_e, u2_t_e] = run_mpc_simulation(A, B, N, Aeq_c, G_c, T_sim, x0, nx, nu, lb_c, ub_c);
y2_t_e = x2_t_e(3, :);

plot_mpc_results(y2_t_e, u2_t_e, N, 'Problem 3e.2: MPC With Varying Input Blocking')











