% signed laplacian (signed ratio cut)

%% basic setting
n = 300; % n = the number of nodes
m = n/2; % m = the size of each community

%% ground truth
Xt = kron(eye(2), ones(m));
Xt(Xt==0)=-1; % Xt = the true cluster matrix

%% parameters of signed SBM
intv = 0.5;
max_param = 35;
max_repeat = 40;

ap = (0 : intv : max_param);
bp = 9; %(0 : intv : max_range);
an = (1 : intv : max_param);
bn = 16; %(1 : intv : max_param);

pp = ap*log(n)/n; % within-community positive
qp = bp*log(n)/n; % across-community positive
pn = an*log(n)/n; % within-community negative
qn = bn*log(n)/n; % across-community negative

%% run
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
                    fprintf('ap = %.1f, bp = %.1f, an = %.1f, bn = %.1f\n', ap(i), bp(j), an(k), bn(l));
                    for repeat = 1 : max_repeat
                        [A_pos, A_neg, ~] = generate_signed_graph(n, pp(i), qp(j), pn(k), qn(l));
                        x_SRC = clustering_signed_graphs_with_signed_laplacian(A_pos, A_neg);
                        dist_SRC = norm(x_SRC*x_SRC'-Xt, 'fro');
                        recovery_tensor(i,j,k,l) = recovery_tensor(i,j,k,l) + (dist_SRC==0);
                    end
                end
            end
        end
    end
end
recovery_rate = recovery_tensor / max_repeat;
elapsedTime = toc;
fprintf("elapsed time is %f\n", elapsedTime);
fprintf("n = %d, max_repeat = %d\n", n, max_repeat);

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
% imshow(apbp, 'InitialMagnification','fit','XData',[0 max_ap],'YData',[0 max_bp]);
%
% colorbar;
% axis on;
% set(gca,'YDir','normal');
% hold on;
%
% cn = 2 - (sqrt(an)-sqrt(bn)).^2;
%
% f =  @(bp,ap)  sqrt(ap) - sqrt(bp) - sqrt(cn);
% fimplicit(f,[0 bp_max 0 ap_max], 'LineWidth', 1.5, 'color', 'r');
% daspect([1 1 1]);
% xlabel('{$\alpha^+$}','Interpreter','latex', 'FontSize', 20, 'LineWidth', 4);
% ylabel('{$\beta^+$}','Interpreter','latex', 'FontSize', 20, 'LineWidth', 4);
% title('GPM');
% hold on
%
% g =  @(bp,ap)  sqrt(bp) - sqrt(ap) - sqrt(cn);
% fimplicit(g,[0 bp_max 0 ap_max], 'LineWidth', 1.5, 'color', 'r');
% daspect([1 1 1]);
% xlabel('{$\alpha^+$}','Interpreter','latex', 'FontSize', 20, 'LineWidth', 4);
% ylabel('{$\beta^+$}','Interpreter','latex', 'FontSize', 20, 'LineWidth', 4);
% title('GPM');

%% alpha- vs beta-
figure();
anbn = zeros(size_an,size_bn);
for k = 1 : size_an
    for l = 1 : size_bn
        anbn(k,l) = recovery_rate(1,1,k,l);
    end
end
imshow(anbn, 'InitialMagnification','fit','XData',[0 max_param],'YData',[0 max_param]);

colorbar;
axis on;
set(gca,'YDir','normal','FontSize',20,'XTick', 5:5:max_param,'YTick', 5:5:max_param);
hold on;

cp = 2 - (sqrt(ap)-sqrt(bp))^2;

f =  @(bn,an)  sqrt(an) - sqrt(bn) - sqrt(cp);
fimplicit(f,[0 max_param 0 max_param], 'LineWidth', 1.5, 'color', 'r');
daspect([1 1 1]);
hold on

g =  @(bn,an)  sqrt(bn) - sqrt(an) - sqrt(cp);
fimplicit(g,[0 max_param 0 max_param], 'LineWidth', 1.5, 'color', 'r');
daspect([1 1 1]);
xlabel('{$\alpha^-$}','Interpreter','latex', 'FontSize', 50, 'LineWidth', 4);
ylabel('{$\beta^-$}','Interpreter','latex', 'FontSize', 50, 'LineWidth', 4);
% title('SPM');

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
set(gca,'YDir','normal','FontSize',20,'XTick', 5:5:max_param,'YTick', 5:5:max_param);
hold on;

f =  @(an,ap)  (sqrt(ap)-sqrt(bp))^2 + (sqrt(an)-sqrt(bn))^2 - 2;
fimplicit(f,[0 max_param 0 max_param], 'LineWidth', 1.5, 'color', 'r');
daspect([1 1 1]);
xlabel('{$\alpha^-$}','Interpreter','latex', 'FontSize', 50, 'LineWidth', 4);
ylabel('{$\alpha^+$}','Interpreter','latex', 'FontSize', 50, 'LineWidth', 4);
% title('SPM');

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