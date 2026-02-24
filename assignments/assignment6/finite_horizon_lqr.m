%% Finite Horizon LQR
clc;

A = [ 1  0.5;
      0  1  ];
b = [ 0.125  0.5]';

Q = diag([2 2]);
R = 2;

[K, P, e] = dlqr(A, b, Q/2, R/2);

fprintf('Optimal Feedback Gain K:\n')
disp(K)

fprintf('Steady-state Riccati Solution P:\n')
disp(P)

fprintf('Closed-loop Eigenvalues:\n')
disp(e)

if all(abs(e) < 1)
    fprintf('Closed-loop system is stable.\n')
else
    fprintf('Closed loop system is unstable.\n')
end
