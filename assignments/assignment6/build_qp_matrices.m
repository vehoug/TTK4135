function [G, Aeq, beq, lb, ub] = build_qp_matrices(A, B, C, x0, Qy, Rt, N, u_min, u_max)
    [nx, nu] = size(B);
    
    Qt = C' * Qy * C;
    Q_big = kron(eye(N), Qt);
    R_big = kron(eye(N), Rt);
    G = blkdiag(Q_big, R_big);
    
    Aeq_x = eye(N*nx) + kron(diag(ones(N-1, 1), -1), -A);
    Aeq_u = kron(eye(N), -B); 
    
    Aeq = [Aeq_x, Aeq_u];
    beq = [A*x0; zeros((N-1)*nx, 1)];
   
    lb = [-Inf(N*nx, 1); u_min * ones(N*nu, 1)];
    ub = [ Inf(N*nx, 1); u_max * ones(N*nu, 1)];
end