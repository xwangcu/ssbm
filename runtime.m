% runtime comparison
warning('off');

%% Basic Setting
n = 5000;       % n = the number of nodes
m = n/2;       % m = the size of each community
max_repeat = 40;

ap = 16;
bp = 9;
an = 9;
bn = 16;
pp = ap*log(n)/n; % within-community positive
qp = bp*log(n)/n; % across-community positive
pn = an*log(n)/n; % within-community negative
qn = bn*log(n)/n; % across-community negative

Xt =  kron(eye(2), ones(m));
Xt(Xt==0)=-1; % Xt = the true cluster matrix

fprintf('n = %d\n', n);
fprintf('ap = %d, bp = %d, an = %d, bn = %d\n', ap, bp, an, bn);
fprintf('max_repeat = %d\n', max_repeat);

%% GPM
ttime_GPM = 0;
for repeat = 1 : max_repeat
    [A_pos, A_neg, ~] = generate_signed_graph(n, pp, qp, pn, qn);
    tic;
    [pp_est, qp_est] = estimate_probability(A_pos);
    ap_est = pp_est / (log(n)/n);
    bp_est = qp_est / (log(n)/n);
    [pn_est, qn_est] = estimate_probability(A_neg);
    an_est = pn_est / (log(n)/n);
    bn_est = qn_est / (log(n)/n);
    xi_est = estimate_weight(ap_est, bp_est, an_est, bn_est);
    A = A_pos - xi_est * A_neg;
    tau = sum(sum(A+(pp-xi_est*pn)*eye(n)))/n^2;
    A = A - tau*ones(n);
    x0 = randn(n,1);
    x0 = x0/norm(x0);
    opts = struct('T', 1e3, 'tol', 1e-4, 'report_interval', 1, 'total_time', 1000); %%% choose the paraneters in GPM
    [x_GPM, iter_GPM] = gpm_ssbm_comp(A, x0, opts);
    time_GPM = toc;
    ttime_GPM = ttime_GPM + time_GPM;
    dist_GPM = norm(x_GPM*x_GPM'-Xt, 'fro');
end
fprintf('GPM time         = %f\n', ttime_GPM);

%% SRC
ttime_SRC = 0;
for repeat = 1 : max_repeat
    [A_pos, A_neg, ~] = generate_signed_graph(n, pp, qp, pn, qn);
    tic;
    x_SRC = clustering_signed_graphs_with_signed_laplacian(A_pos, A_neg);
    time_SRC = toc;
    ttime_SRC = ttime_SRC + time_SRC;
end
fprintf('SRC time         = %f\n', ttime_SRC);

%% SPONGE
tau_pos = 10;
tau_neg = 1;
ttime_SPONGE = 0;
for repeat = 1 : max_repeat
    [A_pos, A_neg, ~] = generate_signed_graph(n, pp, qp, pn, qn);
    tic;
    [x_SPONGE] = clustering_signed_graphs_with_sponge(A_pos, A_neg, tau_pos, tau_neg);
    time_SPONGE = toc;
    ttime_SPONGE = ttime_SPONGE + time_SPONGE;
end
fprintf('SPONGE time      = %f\n', ttime_SPONGE);

%% SPM with p=0;
power = 0;
groundTruth = zeros(n,1);
groundTruth(1:n/2) = 1;
GroundTruthPerLayerCell = {groundTruth, groundTruth};
pinVec = [pp pn];
poutVec = [qp qn];
ttime_SPM1 = 0;
for repeat = 1 : max_repeat
    Wcell = generate_multilayer_graph(2, GroundTruthPerLayerCell, pinVec, poutVec);
    tic;
    x_SPM1 = clustering_signed_graphs_with_power_mean_laplacian(Wcell, power, 2);
    x_SPM1(x_SPM1==1) = 1;
    x_SPM1(x_SPM1==2) = -1;
    time_SPM1 = toc;
    ttime_SPM1 = ttime_SPM1 + time_SPM1;
end
fprintf('SPM (p=0) time   = %f\n', ttime_SPM1);

%% SPM with p=-10;
power = -10;
groundTruth = zeros(n,1);
groundTruth(1:n/2) = 1;
GroundTruthPerLayerCell = {groundTruth, groundTruth};
pinVec = [pp pn];
poutVec = [qp qn];
ttime_SPM2 = 0;
for repeat = 1 : max_repeat
    Wcell = generate_multilayer_graph(2, GroundTruthPerLayerCell, pinVec, poutVec);
    tic;
    x_SPM2 = clustering_signed_graphs_with_power_mean_laplacian(Wcell, power, 2);
    x_SPM2(x_SPM2==1) = 1;
    x_SPM2(x_SPM2==2) = -1;
    time_SPM2 = toc;
    ttime_SPM2 = ttime_SPM2 + time_SPM2;
end
fprintf('SPM (p=-10) time = %f\n', ttime_SPM2);
