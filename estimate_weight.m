function [xi_est] = estimate_weight(ap_est, bp_est,an_est, bn_est)

xi_est = log(bn_est/an_est) / log(ap_est/bp_est);

if isnan(xi_est)==1 || isinf(xi_est)==1
    xi_est = 1;
end

% if (ap_est > bp_est && an_est < bn_est) || (ap_est < bp_est && an_est > bn_est)
%     xi_est = log(bn_est/an_est) / log(ap_est/bp_est);
% else
%     xi_est = 1;
    
end
