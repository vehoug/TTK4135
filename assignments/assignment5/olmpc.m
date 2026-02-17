%% Problem 1: Open Loop Optimal Control of Inverted Pendulum
close all; clc;

set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');

%% System Model
% Inverted pendulum, linearized, discretized with step 0.1s
A = [ 1,  0.1,     0,      0  ;
      0,  0.9818,  0.2673, 0  ;
      0,  0,       1,      0.1;
      0, -0.0455,  3.1182, 1  ];
     
B = [ 0; 0.1818; 0; 0.4545];

C =  [ 1, 0, 0, 0;
       0, 0, 1, 0];
   
x0 = [ 0  0  0.1 0]';

N = 30;

r = 5;

nx = size(A,1); % nx: number of states
nu = size(B,2); % nu: number of controls
ny = size(C,1); % ny: number of controls

% Cost function
Qt = 2 * diag([100, 0, 10, 0]);
Rt = 2 * r;

%% QP formulation

I_N = eye(N);
Q   = kron(I_N, Qt);
R   = kron(I_N, Rt);
G   = blkdiag(Q, R);

% Equality Constraints
Aeq_1 = eye(N*nx);
Aeq_2 = kron(diag(ones(N-1, 1), -1), -A);
Aeq_3 = kron(I_N, -B);

Aeq =  [Aeq_1 + Aeq_2, Aeq_3]; 

beq = ([A*x0; zeros((N-1)*nx, 1)]);

%% Solve as equality-constrained QP via KKT-system (eq. (16.4) in N&W)

% Define and solve KKT system
KKT_matrix = [ G     -Aeq';
               Aeq    zeros(N*nx)];
KKT_vector = [ zeros(N*(nx+nu), 1);  beq];

KKT_solution = KKT_matrix \ KKT_vector;

% Extract variables
z1 = KKT_solution(1:N*(nx+nu));
plotinvpend(z1, N, x0, 'Open Loop OC: KKT-system Solution');

%% Solve equality-constrained QP with quadprog

opt = optimset('Display','iter-detailed','Algorithm', 'active-set');
[z2, ~, ~, ~, ~]  = quadprog(G, [], [], [], Aeq, beq, ...
                            [], [], zeros(1, N*(nx+nu)), opt);
plotinvpend(z2, N, x0, 'Open Loop OC: Eq. Constrained QP');

%% Solve using quadprog with inequality constraint

% Define inequality constraints
x_lb = -Inf(N*nx, 1);
x_ub =  Inf(N*nx, 1);
u_lb = -ones(N*nu, 1);
u_ub =  ones(N*nu, 1);
lb   =  [x_lb; u_lb];
ub   =  [x_ub; u_ub];

opt = optimset('Display','iter-detailed','Algorithm', 'active-set');
[z3, ~, ~, ~, ~] = quadprog (G, [], [], [], Aeq, beq,...
                            lb, ub, zeros(1, N*(nx+nu)), opt); 

plotinvpend(z3, N, x0, 'Open Loop OC: Ineq. Constrained QP');

x_vec = z3(1 : N*nx);
u_vec = z3(N*nx + 1 : end);

x_ol = reshape(x_vec, nx, N);
x_ol = [x0, x_ol];
u_ol = reshape(u_vec, nu, N);

%% Problem 2b: MPC (Horizon Length Comparisons)

T_sim = 30;

% MPC with N=30
[G30, Aeq30, lb30, ub30] = build_qp_matrices(A, B, Qt, Rt, 30, nx, nu);
[x_mpc30, u_mpc30] = run_mpc_simulation(A, B, 30, Aeq30, G30, T_sim, x0, nx, nu, lb30, ub30);

figure(2)
compare_trajectories(x_ol, u_ol, x_mpc30, u_mpc30, 'Open Loop OC', 'MPC ($N=30$)', 'Open Loop OC vs. MPC w. Long Horizon')

% MPC with N=10
[G10, Aeq10, lb10, ub10] = build_qp_matrices(A, B, Qt, Rt, 10, nx, nu);

[x_mpc10, u_mpc10] = run_mpc_simulation(A, B, 10, Aeq10, G10, T_sim, x0, nx, nu, lb10, ub10);

figure(3)
compare_trajectories(x_ol, u_ol, x_mpc10, u_mpc10, 'Open Loop OC', 'MPC ($N=10$)', 'Open Loop OC vs. MPC w. Short Horizon')

%% Problem 2c: MPC (Plant and Model Mismatch)

B_plant = 0.98 * B;
[x_mism, u_mism] = run_mpc_simulation(A, B_plant, 30, Aeq30, G30, T_sim, x0, nx, nu, lb30, ub30);
compare_trajectories(x_mpc30, u_mpc30, x_mism, u_mism, 'MPC (Perfect Model)', 'MPC (Model/Plant Mismatch)', 'Model Mismatch Effect on MPC ($N=30$)')

%% Helper Functions

function [G, Aeq, lb, ub] = build_qp_matrices(A, B, Qt, Rt, N, nx, nu)
    I_N = eye(N);
    Q_big = kron(I_N, Qt);
    R_big = kron(I_N, Rt);
    G = blkdiag(Q_big, R_big);

    Aeq_x = eye(N*nx) + kron(diag(ones(N-1, 1), -1), -A);
    
    Aeq_u = kron(I_N, -B);
    
    Aeq = [Aeq_x, Aeq_u];
   
    x_lb = -Inf(N*nx, 1);
    x_ub =  Inf(N*nx, 1);
    u_lb = -ones(N*nu, 1);
    u_ub =  ones(N*nu, 1);
    
    lb = [x_lb; u_lb];
    ub = [x_ub; u_ub];
end

function plotinvpend(z, N, x0, title)
    nx = 4;
    nu = 1;
    
    pos = [x0(1); z(1:nx:N*nx-3)];
    phi = [x0(3); z(3:nx:N*nx-1)];
    u = z(N*nx+1:N*nx+N*nu);
    
    t = 1:N;
    
    figure;
    sgtitle(title)
    subplot(3,1,1);
    plot([0,t],pos,'-ro', 'LineWidth', 1.5);
    grid('on');
    ylabel('$\chi_t$', 'FontSize', 14)
    subplot(3,1,2);
    plot([0,t],phi,'-bo', 'LineWidth', 1.5);
    grid('on');
    ylabel('$\phi_t$', 'FontSize', 14)
    subplot(3,1,3);
    stairs(t-1,u,'-mo', 'LineWidth', 1.5);
    grid('on');
    xlabel('$t$', 'FontSize', 14);
    ylabel('$u_t$', 'FontSize', 14);
end

function compare_trajectories(x1, u1, x2, u2, name1, name2, main_title)
    t1 = 0:size(x1, 2)-1;
    t2 = 0:size(x2, 2)-1;
    tu1 = 0:size(u1, 2)-1;
    tu2 = 0:size(u2, 2)-1;

    sgtitle(main_title, 'Interpreter', 'latex');
    
    subplot(3,1,1);
    plot(t1, x1(1,:), '-bo', 'LineWidth', 1.5); hold on;
    plot(t2, x2(1,:), ':ro', 'LineWidth', 1.5);
    ylabel('$\chi_t$', 'FontSize', 14); grid on;
    legend(name1, name2, 'Location', 'best');
    
    subplot(3,1,2);
    plot(t1, x1(3,:), '-bo', 'LineWidth', 1.5); hold on;
    plot(t2, x2(3,:), ':ro', 'LineWidth', 1.5);
    ylabel('$\phi_t$', 'FontSize', 14); grid on;
    
    subplot(3,1,3);
    stairs(tu1, u1(1,:), '-bo', 'LineWidth', 1.5); hold on;
    stairs(tu2, u2(1,:), ':ro', 'LineWidth', 1.5);
    ylabel('$u_t$', 'FontSize', 14); xlabel('$t$', 'FontSize', 14); grid on;
end

function [x_hist, u_hist] = run_mpc_simulation(A, B, N, Aeq, G, T_sim, x0, nx, nu, lb, ub)
    options = optimoptions('quadprog','Display','off');
    
    x_t = x0;
    x_hist = zeros(nx, T_sim + 1);
    u_hist = zeros(nu, T_sim);
    x_hist(:, 1) = x0;
    
    for t = 1:T_sim
        beq = [A * x_t; zeros((N-1)*nx, 1)];
        
        [z, ~, ~] = quadprog(G, [], [], [], Aeq, beq, lb, ub, [], options);
        
        u_start_idx = N*nx + 1;
        u_apply = z(u_start_idx : u_start_idx + nu - 1);
        
        x_next = A * x_t + B * u_apply;
        
        u_hist(:, t) = u_apply;
        x_hist(:, t+1) = x_next;
        
        x_t = x_next;
    end
end