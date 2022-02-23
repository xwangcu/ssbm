function [cluster_ind] = clustering_signed_graphs_with_signed_laplacian(A_pos, A_neg)

% compute the signed laplacian matrix
A = A_pos - A_neg;
D_signed = diag(sum(abs(A),2));
L_signed = D_signed - A;

% compute the eigenvector corresponding to the smallest eigenvalue
[EigVec, EigVal] = eig(L_signed);
[eigVal, eigVal_ind] = sort(diag(EigVal));

% k-means clustering
eps = 1e-10;
isPSD = all(eigVal>=-eps);
if isPSD == 1
    cluster_ind = sign(EigVec(:,eigVal_ind(1)));
end