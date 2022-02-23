function [cluster_id] = clustering_signed_graphs_with_sponge(A_pos, A_neg, tau_pos, tau_neg)

% compute the signed laplacian matrix
A = A_pos - A_neg;
D_pos = diag(sum(abs(A_pos),2));
L_pos = D_pos - A;
D_neg = diag(sum(abs(A_neg),2));
L_neg = D_neg - A;

% compute the bottom 2 generalized eigenvectors
L1 = L_pos + tau_neg*D_neg; % L^+ + tau^-*D^-
L2 = L_neg + tau_pos*D_pos; % L^- + tau^+*D^+
[W,D] = eig(L1,L2);
[~, eigInd] = sort(diag(D));
W2 = W(:,[eigInd(1),eigInd(2)]);

[V_L2,D_L2] = eig(L2);
D_L2(D_L2>0) = sqrt(D_L2(D_L2>0));
L2_sr = V_L2 * D_L2 * V_L2';
V2 = L2_sr * W2;
% v2 = real(sqrtm(L2))*w2;

% clustering
cluster_id = kmeans(V2,2);
cluster_id(cluster_id==1) = 1;
cluster_id(cluster_id==2) = -1;
