function [G, Aeq, beq, lb, ub, T] = build_qp_blocked(A, B, C, x0, Qy, Rt, block_lengths)
    nx = size(A, 1);
    nu = size(B, 2);
    N = sum(block_lengths);
    nb = length(block_lengths);

    T_cell = cell(1, nb);
    for i = 1:nb
        T_cell{i} = ones(block_lengths(i), 1);
    end
    T = blkdiag(T_cell{:});
    
    Qt = C' * Qy * C;
    Q_big = kron(eye(N), Qt);
    
    R_big = kron(eye(N), Rt);
    G_v = T' * R_big * T; 
    G = blkdiag(Q_big, G_v);

    Aeq_x = eye(N*nx) + kron(diag(ones(N-1, 1), -1), -A);
    
    Aeq_u = kron(T, -B); 
    
    Aeq = [Aeq_x, Aeq_u];
    beq = [A*x0; zeros((N-1)*nx, 1)];

    lb = [-Inf(N*nx, 1); -ones(nb*nu, 1)];
    ub = [ Inf(N*nx, 1);  ones(nb*nu, 1)];
end
