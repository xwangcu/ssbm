function [p_est, q_est] = estimate_probability(A)

%% extract graph information
n = size(A,1);
m = n/2;

%% compute the number of edges and triangles (up to a scalar factor)
E = sum(sum(A));
T = trace(A*A*A);

%% estimate connecting probabilities
cubic_para = [m^2-3*m+2+3*(m-1)^3/m, -3*(m-1)^2*E/m^2, 3*(m-1)*E^2/(4*m^3), -T/n];
sol = roots(cubic_para); 
p_est = sol(imag(sol)==0);
if size(p_est,1) > 1
    p_est = p_est(1);
end
q_est = (E-(n^2/2-n)*p_est) / (n^2/2);

