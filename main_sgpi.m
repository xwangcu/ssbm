% signed SBM with estimated parameters alpha+ & beta+ & alpha- & beta-

%% Basic Setting
n = 300;       % n = the number of nodes
m = n/2;       % m = the size of each community
Xt =  kron(eye(2), ones(m));
Xt(Xt==0)=-1;                           %%% Xt = the true cluster matrix
xt = [ones(m,1); -ones(m,1)];           %%% xt = the true cluster vector

%% Parameters of Signed SBM
intv = 0.5;
max_param = 35;
max_repeat = 40;

ap = (1 : intv : max_param);
bp = 9; %(1 : intv : max_param);
an = (1 : intv : max_param);
bn = 16; %(1 : intv : max_param);

pp = ap*log(n)/n; % within-community positive
qp = bp*log(n)/n; % across-community positive
pn = an*log(n)/n; % within-community negative
qn = bn*log(n)/n; % across-community negative

%% Run Algorithm
ttime_GPM = deal(0);
size_ap = size(ap,2);
size_bp = size(bp,2);
size_an = size(an,2);
size_bn = size(bn,2);
recovery_tensor = zeros(size_ap,size_bp,size_an,size_bn);
tic
for i = 1 : size_ap
    for j = 1 : size_bp
        for k = 1 : size_an
            for l = 1 : size_bn
                if 1 || (ap(i) ~= bp(j) && an(k) ~= bn(l))
                    fprintf('ap = %.2f, bp = %.2f, an = %.2f, bn = %.2f \n', ap(i), bp(j), an(k), bn(l));

                    % initialization
                    x0 = randn(n,1);
                    x0 = x0/norm(x0);
                    xi = log(bn(l)/an(k)) / log(ap(i)/bp(j));
                    
                    % repeat
                    for repeat = 1 : max_repeat
                        % generate random graph within community 1
                        Ans11 = rand(m);
                        Al11 = tril(Ans11,-1);
                        As11 = Al11 + Al11';
                        A11_pos = double(As11<pp(i)) - diag(ones(m,1));
                        A11_neg = double(As11>=pp(i) & As11<=pp(i)+pn(k));
                        
                        % generate random graph within community 2
                        Ans22 = rand(m);
                        Al22 = tril(Ans22,-1);
                        As22 = Al22 + Al22';
                        A22_pos = double(As22<pp(i)) - diag(ones(m,1));
                        A22_neg = double(As22>=pp(i) & As22<=pp(i)+pn(k));
                        
                        % generate random graph aross communities
                        As12 = rand(m);
                        A12_pos = double(As12<qp(j));
                        A12_neg = double(As12>=qp(j) & As12<=qp(j)+qn(l));
                        
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
                        A = ([A11, A12; A12', A22]);
                        tau = sum(sum(A))/n^2;
                        A = A - tau*ones(n);
                        
                        % run
                        opts = struct('T', 1e3, 'tol', 1e-4, 'report_interval', 1, 'total_time', 1000); %%% choose the paraneters in GPM
                        tic;
                        [x_GPM, iter_GPM, val_collector_GPM, itergap_GPM] = gpm_ssbm(A, x0, xt, opts);
                        time_GPM = toc;
                        ttime_GPM = ttime_GPM + time_GPM;
                        dist_GPM = norm(x_GPM*x_GPM'-Xt, 'fro');
                        recovery_tensor(i,j,k,l) = recovery_tensor(i,j,k,l) + (dist_GPM==0);
                    end
                end
            end
        end
    end
end
recovery_rate = recovery_tensor / max_repeat;
elapsedTime = toc;
fprintf("elapsed time is %f\n", elapsedTime);

%% alpha+ vs beta+
% figure();
% apbp = zeros(size_ap,size_bp);
% ind_an = find(an==3);
% ind_bn = find(bn==4);
% for i = 1 : size_ap
%     for j = 1 : size_bp
%         apbp(i,j) = recovery_rate(i,j,ind_an,ind_bn);
%     end
% end
% imshow(apbp, 'InitialMagnification','fit','XData',[0 max_param],'YData',[0 max_param]);
% 
% colorbar;
% axis on;
% set(gca,'YDir','normal');
% hold on;
% 
% cn = 2 - (sqrt(an)-sqrt(bn)).^2;
% 
% f =  @(bp,ap)  sqrt(ap) - sqrt(bp) - sqrt(cn);
% fimplicit(f,[0 max_param 0 max_param], 'LineWidth', 1.5, 'color', 'r');
% daspect([1 1 1]);
% xlabel('{$\alpha^+$}','Interpreter','latex', 'FontSize', 20, 'LineWidth', 4);
% ylabel('{$\beta^+$}','Interpreter','latex', 'FontSize', 20, 'LineWidth', 4);
% title('GPM');
% hold on
% 
% g =  @(bp,ap)  sqrt(bp) - sqrt(ap) - sqrt(cn);
% fimplicit(g,[0 max_param 0 max_param], 'LineWidth', 1.5, 'color', 'r');
% daspect([1 1 1]);
% xlabel('{$\alpha^+$}','Interpreter','latex', 'FontSize', 20, 'LineWidth', 4);
% ylabel('{$\beta^+$}','Interpreter','latex', 'FontSize', 20, 'LineWidth', 4);
% title('GPM');

%% alpha- vs beta-
% figure();
% anbn = zeros(size_an,size_bn);
% for k = 1 : size_an
%     for l = 1 : size_bn
%         anbn(k,l) = recovery_rate(1,1,k,l);
%     end
% end
% anbn = anbn(1:69, 1:69);
% imshow(anbn, 'InitialMagnification','fit','XData',[1 35],'YData',[1 35]);
% colorbar;
% axis on;
% set(gca,'YDir','normal','FontSize',20,'XTick', 5:5:35);
% hold on;
% 
% cp = 2 - (sqrt(ap)-sqrt(bp))^2;
% f =  @(bn,an)  sqrt(an) - sqrt(bn) - sqrt(cp);
% fimplicit(f,[0 35 0 35], 'LineWidth', 1.5, 'color', 'r');
% daspect([1 1 1]);
% hold on
% g =  @(bn,an)  sqrt(bn) - sqrt(an) - sqrt(cp);
% fimplicit(g,[0 35 0 35], 'LineWidth', 1.5, 'color', 'r');
% daspect([1 1 1]);
% xlabel('{$\alpha^-$}','Interpreter','latex', 'FontSize', 50, 'LineWidth', 4);
% ylabel('{$\beta^-$}','Interpreter','latex', 'FontSize', 50, 'LineWidth', 4);
% % title('SGPM');

%% alpha+ vs alpha-
figure();
apan = zeros(size_ap,size_an);
for i = 1 : size_ap
    for k = 1 : size_an
        apan(i,k) = recovery_rate(i,1,k,1);
    end
end
imshow(apan, 'InitialMagnification','fit','XData',[1 max_param],'YData',[1 max_param]);
colorbar;
axis on;
set(gca,'YDir','normal','FontSize',20,'XTick', 5:5:35);
hold on;

f =  @(an,ap)  (sqrt(ap)-sqrt(bp))^2 + (sqrt(an)-sqrt(bn))^2 - 2;
fimplicit(f,[0 max_param 0 max_param], 'LineWidth', 1.5, 'color', 'r');
daspect([1 1 1]);
xlabel('{$\alpha^-$}','Interpreter','latex', 'FontSize', 50, 'LineWidth', 4);
ylabel('{$\alpha^+$}','Interpreter','latex', 'FontSize', 50, 'LineWidth', 4);
title('GPM');

%% beta+ vs beta-
% figure();
% bpbn = zeros(size_bp,size_bn);
% for j = 1 : size_bp
%     for l = 1 : size_bn
%         bpbn(j,l) = recovery_rate(1,j,1,l);
%     end
% end
% imshow(bpbn, 'InitialMagnification','fit','XData',[1 max_param],'YData',[1 max_param]);
% colorbar;
% axis on;
% set(gca,'YDir','normal');
% hold on;
% 
% f =  @(bn,bp)  (sqrt(ap)-sqrt(bp))^2 + (sqrt(an)-sqrt(bn))^2 - 2;
% fimplicit(f,[0 max_param 0 max_param], 'LineWidth', 1.5, 'color', 'r');
% daspect([1 1 1]);
% xlabel('{$\beta^-$}','Interpreter','latex', 'FontSize', 20, 'LineWidth', 4);
% ylabel('{$\beta^+$}','Interpreter','latex', 'FontSize', 20, 'LineWidth', 4);
% title('GPM');
