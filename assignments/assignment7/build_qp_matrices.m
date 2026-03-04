function [G, Aeq, beq, lb, ub] = build_qp_matrices(A, B, x0, Q, R, N, u_min, u_max)
    [nx, nu] = size(B);
    
    Q_big = kron(eye(N), Q);
    R_big = kron(eye(N), R);
    G = blkdiag(Q_big, R_big);
    
    Aeq_x = eye(N*nx) + kron(diag(ones(N-1, 1), -1), -A);
    Aeq_u = kron(eye(N), -B); 
    
    Aeq = [Aeq_x, Aeq_u];
    
    beq = [A*x0; zeros((N-1)*nx, 1)];
   
    lb = [-Inf(N*nx, 1); u_min * ones(N*nu, 1)];
    ub = [ Inf(N*nx, 1); u_max * ones(N*nu, 1)];
end

