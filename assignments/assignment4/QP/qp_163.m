function [x, fval, iterates, W, lambda, foo] = qp_163(G,c,A,b,x0,W0,outp)
%QP_163 Quadratic programming based on algorithm 16.3 in Nocedal&Wright
%   [X,FVAL,ITERATES, W, LAMBDA, FOO] = qp_163(G,c,A,b,x0,W0,outp) attempts
%   to solve the quadratic programming problem:
%        
%            min 0.5*x'*G*x + c'*x    subject to:   A*x >= b 
%             x
%   
%   The starting point is set by x0, and this point has to be feasible, 
%   with an associated (subset of) active constraints W0 (a vector the size 
%   of number of rows of A, with 1 for the active constraints, and 0 for 
%   the rest). 
% 
%   If outp is passed as 'report', the function prints a detailed report at
%   every iteration. If report is passed as [], no report is printed. At
%   termination, X contains the optimal point, FVAL the optimal objective
%   function, and ITERATES is a matrix where each column i contains x at 
%   iteration i.
%
%   The output W is the active set corresponding to x, and LAMBDA is the
%   Lagrange multipliers for the active constraints.
%
%   The output FOO is an evaluation of the first order optimality (KKT
%   stationarity) condition, and should ideally be zero upon convergence.
%   
%   NOTE: This function is written for educational purposes only. The
%         algorithm is neither quick nor robust, and should only be used on
%         small problems with feasible initial conditions. 
%         The implementation is meant to illustrate Algorithm 16.3 in 
%         Nocedal and Wright: Numerical Optimization, 2nd ed., Springer, 
%         2006. There are probably some errors in the
%         implementation, please report any to lecturer of TTK4135.
% 
%   Version 0.1
%   Lars Imsland, 07.11.2024

% Note: We do not do any sanity check on the inputs

n = size(G,1);
m = size(A,1);

I = eye(n); % Identity matrix

% A matrix containing x for all iterations
iterates(:,1) = x0;

Wkp = W0; 

invG = inv(G);

k=0; % iteration counter
opt = 0; % Optimum not found
blockind = 0; alpha = 1; % Define beforehand for reporting

while not(opt),
    k = k+1;
    xk = iterates(:,k);
    Wk = Wkp;
    Ak = A(find(Wk),:);
    bk = b(find(Wk));
    gk = G*xk + c;
    hk = Ak*xk - bk;
    nA = size(Ak,1); % Number of active constraints in working set

    % Solve (16.39) using direct solution of KKT-system (16.5)
    % Note that this can be done more efficiently, see 16.2 and 16.3 in
    % Nocedal&Wright
    
    z = linsolve([G, Ak';Ak, zeros(nA)],[gk; hk]);
    pk = -z(1:n);
    lambda = z(n+1:end);

    if norm(pk) < 1e-8, % if pk == 0
        % we have already calculated lambdahat (16.42) while solving
        % KKT-system above
        if all(lambda >= 0), % if all elements of lambda are positive
            % We have found a point fulfilling KKT conditions
            x = xk;
            opt = 1;
        else
            % find index of smallest lagrange multiplier
            [minlam,jl] = min(lambda); % jl is index in lambda, 
                                      % need to convert to index in Wk
            tmp = find(Wk);
            j = tmp(jl);
            Wkp(j) = 0; % Remove active constraint
            iterates(:,k+1) = xk;
        end
    else % pk is not zero
        [alpha,blockind] = find_alpha(A,b,xk,Wk,pk);
        iterates(:,k+1) = xk + alpha*pk;
        if blockind>0,
            Wkp(blockind) = 1;
        % else Wk stays unchanged
        end  
    end

    if strcmpi(outp, 'report')
        report(k, xk, Wk, pk, lambda, blockind, alpha)
    end
end

fval = .5*x'*G*x + c'*x;
W = Wkp;

% check first order optimality (stationarity KKT)
foo = norm(G*x + c - Ak'*lambda);

if opt == 1,
    finalreport(x, fval);
end
end


function [alpha, blockind] = find_alpha(A,b,xk,Wk,pk)
alpha = 1;
blockind = 0;
for i = 1:size(A,1),
    if (Wk(i) == 0) & (A(i,:)*pk < 0),
        alpha_tmp = (b(i)-A(i,:)*xk)/(A(i,:)*pk);
        if alpha_tmp < alpha,
            alpha = alpha_tmp;
            blockind = i;
        end
    end
end
end

function [] = report(iteration, xk, Wk, pk, lambda, blockind, alpha)
    line = repmat('-',1,80);
    fprintf('\n%s\n',line);
    fprintf('    Iteration number: %i\n', iteration-1);
    fprintf(['    Current iterate xk           [' sprintf('%9.4f, ', xk) '\b\b]''\n' ]);
    fprintf(['    Current working set Wk       {' sprintf('%i, ', Wk) '\b\b}\n' ]);
    fprintf(['    pk                           [' sprintf('%9.4f, ', pk) '\b\b]''\n' ]);
    fprintf(['    lambda (for working set)     [' sprintf('%9.4f, ', lambda) '\b\b]''\n' ]);
    if norm(pk) > 1e-8
        fprintf(['    alpha                         ' sprintf('%i, ', alpha) '\b\b\n' ]);
        fprintf(['    blocking index (0 if none)   {' sprintf('%i, ', blockind) '\b\b}\n' ]);
    end
end

function [] = finalreport(x, fval)
    line = repmat('-',1,80);
    fprintf('\n');
    fprintf('    OPTIMAL POINT FOUND\n');
    fprintf(['    x^* =    [' sprintf('%9.4f, ', x) '\b\b]''\n' ]);
    fprintf( '    0.5 x''G x + c''x^* = %14.8f\n', fval);
    fprintf('%s\n\n',line);
end

