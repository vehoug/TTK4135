function [x_hist, u_hist] = run_mpc_simulation(A, B, N, Aeq, G, T_sim, x0, nx, nu, lb, ub)
    options = optimoptions('quadprog', 'Display', 'off');
    
    x_t = x0;
    x_hist = zeros(nx, T_sim + 1);
    u_hist = zeros(nu, T_sim);
    x_hist(:, 1) = x0;
    
    for t = 1:T_sim
        beq = [A * x_t; zeros((N-1)*nx, 1)];
        
        [z, ~, exitflag] = quadprog(G, [], [], [], Aeq, beq, lb, ub, [], options);
        
        if exitflag ~= 1
            warning('Solver failed at time step %d with exitflag %d', t, exitflag);
        end
        
        u_start_idx = N*nx + 1;
        u_apply = z(u_start_idx : u_start_idx + nu - 1);
        
        x_next = A * x_t + B * u_apply;
        
        u_hist(:, t) = u_apply;
        x_hist(:, t+1) = x_next;
        
        x_t = x_next;
    end
end

