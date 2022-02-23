% runtime comparison
warning('off');

%% Load Real Graph
realGraphName = 'congress';
Adj = read_real_graph_data(realGraphName);
colorbar
A_pos = double(Adj>0);
A_neg = double(Adj<0);
n = size(Adj,1); % n = the number of nodes
m = n/2;       % m = the number of nodes in each community
fprintf('n = %d\n', n);

%% GPM
tic;
[pp_est, qp_est] = estimate_probability(A_pos);
ap_est = pp_est / (log(n)/n);
bp_est = qp_est / (log(n)/n);
[pn_est, qn_est] = estimate_probability(A_neg);
an_est = pn_est / (log(n)/n);
bn_est = qn_est / (log(n)/n);
xi_est = estimate_weight(ap_est, bp_est, an_est, bn_est);
A = A_pos - xi_est * A_neg;
tau = sum(sum(A+(pp_est-xi_est*pn_est)*eye(n)))/n^2;
A = A - tau*ones(n);
x0 = randn(n,1);
x0 = x0/norm(x0);
opts = struct('T', 1e3, 'tol', 1e-4, 'report_interval', 1, 'total_time', 1000);
[x_GPM, iter_GPM] = gpm_ssbm_comp(A, x0, opts);
time_GPM = toc;
fprintf('GPM time         = %f\n', time_GPM);

%% SRC
tic;
x_SRC = clustering_signed_graphs_with_signed_laplacian(A_pos, A_neg);
time_SRC = toc;
fprintf('SRC time         = %f\n', time_SRC);

%% SPONGE
tau_pos = 10;
tau_neg = 1;
tic;
[x_SPONGE] = clustering_signed_graphs_with_sponge(A_pos, A_neg, tau_pos, tau_neg);
time_SPONGE = toc;
fprintf('SPONGE time      = %f\n', time_SPONGE);

%% SPM with p=-10
power = -10;
Wcell = {A_pos, A_neg};
tic;
x_SPM2 = clustering_signed_graphs_with_power_mean_laplacian(Wcell, power, 2);
x_SPM2(x_SPM2==1) = 1;
x_SPM2(x_SPM2==2) = -1;
time_SPM2 = toc;
fprintf('SPM (p=-10) time = %f\n', time_SPM2);

%% Plot the Sorted Adjacency Matrices
ind1 = find(x_GPM==1);
ind2 = find(x_GPM==-1);
Row1 = Adj(ind1,:);
Row2 = Adj(ind2,:);
C11_GPM = Row1(:,ind1);
C12_GPM = Row1(:,ind2);
C21_GPM = Row2(:,ind1);
C22_GPM = Row2(:,ind2);
C_GPM = [C11_GPM,C12_GPM;C21_GPM,C22_GPM];

ind1 = find(x_SRC==1);
ind2 = find(x_SRC==-1);
Row1 = Adj(ind1,:);
Row2 = Adj(ind2,:);
C11_SRC = Row1(:,ind1);
C12_SRC = Row1(:,ind2);
C21_SRC = Row2(:,ind1);
C22_SRC = Row2(:,ind2);
C_SRC = [C11_SRC,C12_SRC;C21_SRC,C22_SRC];

ind1 = find(x_SPONGE==1);
ind2 = find(x_SPONGE==-1);
Row1 = Adj(ind1,:);
Row2 = Adj(ind2,:);
C11_SPONGE = Row1(:,ind1);
C12_SPONGE = Row1(:,ind2);
C21_SPONGE = Row2(:,ind1);
C22_SPONGE = Row2(:,ind2);
C_SPONGE = [C11_SPONGE,C12_SPONGE;C21_SPONGE,C22_SPONGE];

ind1 = find(x_SPM2==1);
ind2 = find(x_SPM2==-1);
Row1 = Adj(ind1,:);
Row2 = Adj(ind2,:);
C11_SPM2 = Row1(:,ind1);
C12_SPM2 = Row1(:,ind2);
C21_SPM2 = Row2(:,ind1);
C22_SPM2 = Row2(:,ind2);
C_SPM2 = [C11_SPM2,C12_SPM2;C21_SPM2,C22_SPM2];

colormap parula
subplot(1,5,1)
imagesc(Adj)
title('Original Adjacency')
subplot(1,5,2)
imagesc(C_GPM)
title('SGPM')
subplot(1,5,3)
imagesc(C_SRC)
title('SRC')
subplot(1,5,4)
imagesc(C_SPONGE)
title('SPONGE')
subplot(1,5,5)
imagesc(C_SPM2)
title('SPM (p=-10)')

%% Function Values
fval_GPM = x_GPM'*A*x_GPM;
fval_SRC = x_SRC'*A*x_SRC;
fval_SPONGE = x_SPONGE'*A*x_SPONGE;
fval_SPM2 = x_SPM2'*A*x_SPM2;

fprintf('fval_GPM         = %d\n', fval_GPM);
fprintf('fval_SRC         = %d\n', fval_SRC);
fprintf('fval_SPONGE      = %d\n', fval_SPONGE);
fprintf('fval_SPM2        = %d\n', fval_SPM2);

%% Metrics
% xxt_GPM = x_GPM * x_GPM';
% Nin_pos_GPM = sum(sum(xxt_GPM==A_pos))/2;
% Nout_neg_GPM = sum(sum(xxt_GPM==-A_neg))/2;
% 
% xxt_SRC = x_SRC * x_SRC';
% Nin_pos_SRC = sum(sum(xxt_SRC==A_pos))/2;
% Nout_neg_SRC = sum(sum(xxt_SRC==-A_neg))/2;
% 
% xxt_SPONGE = x_SPONGE * x_SPONGE';
% Nin_pos_SPONGE = sum(sum(xxt_SPONGE==A_pos))/2;
% Nout_neg_SPONGE = sum(sum(xxt_SPONGE==-A_neg))/2;
% 
% xxt_SPM2 = x_SPM2 * x_SPM2';
% Nin_pos_SPM2 = sum(sum(xxt_SPM2==A_pos))/2;
% Nout_neg_SPM2 = sum(sum(xxt_SPM2==-A_neg))/2;
% 
% fprintf('N+in:\n');
% fprintf('GPM         = %d\n', Nin_pos_GPM);
% fprintf('SRC         = %d\n', Nin_pos_SRC);
% fprintf('SPONGE      = %d\n', Nin_pos_SPONGE);
% fprintf('SPM (p=-10) = %d\n', Nin_pos_SPM2);
% 
% fprintf('N-out:\n');
% fprintf('GPM         = %d\n', Nout_neg_GPM);
% fprintf('SRC         = %d\n', Nout_neg_SRC);
% fprintf('SPONGE      = %d\n', Nout_neg_SPONGE);
% fprintf('SPM (p=-10) = %d\n', Nout_neg_SPM2);
