function [x, iter] = gpm_ssbm_comp(A, x, opts)

%%  GPM solves regularized MLE
% --- INPUT ---
% A:    adjacency matrix (a sparse 0-1 matrix).
% x:    starting point
% xt:   ground truth
% opts.tol:    desired tolerance for suboptimality
% opt.T:    maximal number of iterations

% --- OUTPUT ---
% x: returned clusters by GPM
% iter: terminal number of iteration

%% Parameter setting
maxiter = opts.T;
tol = opts.tol;
n = size(A, 1);
stage1 = floor(0.5*log(n)/log(log(n)));
Ax = A*x; %%% matrix-vector product

for iter = 1 : maxiter
    xold = x;
    
    %% Check fixed-point condition
    x1 = x + Ax;  
    x1(x1>=0) = 1; 
    x1(x1<0) = -1;
    dist = norm(x - x1);
    
    %% update of GPM + PGD
    if iter <= stage1
        % power iteration
        x = Ax;
        x = x/norm(x)*sqrt(n);
    else
        % GPM iteration
        x = Ax;
        x(x>=0) = 1;
        x(x<0) = -1;
    end
    Ax = A*x;

    %%  stopping criterion
    if dist <= tol && (norm(x-xold) == 0 || norm(x+xold) == 0)
        break;
    end    
end
