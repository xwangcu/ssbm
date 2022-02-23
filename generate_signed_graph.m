function [A_pos, A_neg, Atilde] = generate_signed_graph(n, pp, qp, pn, qn)

m = n/2;

% generate random graph within community 1
Ans11 = rand(m);
Al11 = tril(Ans11,-1);
As11 = Al11 + Al11';
A11_pos = double(As11<pp) - diag(ones(m,1));
A11_neg = double(As11>=pp & As11<=pp+pn);

% generate random graph within community 2
Ans22 = rand(m);
Al22 = tril(Ans22,-1);
As22 = Al22 + Al22';
A22_pos = double(As22<pp) - diag(ones(m,1));
A22_neg = double(As22>=pp & As22<=pp+pn);

% generate random graph aross communities
As12 = rand(m);
A12_pos = double(As12<qp);
A12_neg = double(As12>=qp & As12<=qp+qn);

% adjacency matrices of positive & negative graphs
A_pos = ([A11_pos, A12_pos; A12_pos', A22_pos]);
A_neg = ([A11_neg, A12_neg; A12_neg', A22_neg]);

% estimate SSBM parameters
[pp_est, qp_est] = estimate_probability(A_pos);
ap_est = pp_est / (log(n)/n);
bp_est = qp_est / (log(n)/n);
[pn_est, qn_est] = estimate_probability(A_neg);
an_est = pn_est / (log(n)/n);
bn_est = qn_est / (log(n)/n);
xi_est = estimate_weight(ap_est, bp_est, an_est, bn_est);

A11 = A11_pos - xi_est * A11_neg;
A22 = A22_pos - xi_est * A22_neg;
A12 = A12_pos - xi_est * A12_neg;
Atilde = ([A11, A12; A12', A22]);
