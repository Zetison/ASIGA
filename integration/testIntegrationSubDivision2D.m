close all
% 
% rule = 'Freud';
rule = 'Legendre';
N = 200;
n_qpSub = 5;
I_gauss = zeros(N,1);
N_arr = 1:N;
p = 10;
epsilon = 2/(2*p);
% f = @(x,y) 1./(1+(x.^2+y.^2)/epsilon.^2);
f = @(x,y) 1./sqrt((x-0.5).^2+(y-(1+epsilon)).^2);
x = linspace(-1,1,1000);
y = linspace(-1,1,1000);
[X,Y] = meshgrid(x,y);
surf(X,Y,f(X,Y),'EdgeColor','none','LineStyle','none')
return
for n = N_arr
    [Q1D,W1D] = getQuadFromFile(n);
    Q_xi = repmat(Q1D,n,1);
    Q_eta = repmat(Q1D.',n,1);
    Q_eta = Q_eta(:);
    W2D_1 = W1D*W1D.';
    W2D_1 = W2D_1(:);

    I_gauss(n) = sum(f(Q_xi,Q_eta).*W2D_1);
end
I = integral2(f,mp('-1'),mp('1'),mp('-1'),mp('1'),'AbsTol', 1e-20,'RelTol', 1e-20);
[Q1D,W1D] = getQuadFromFile(n_qpSub);
W = W1D*W1D.';
W = W(:);
Q_xi = repmat(Q1D,n_qpSub,1);
Q_eta = repmat(Q1D.',n_qpSub,1);
Q_eta = Q_eta(:);
n_div_arr = 1:round(N_arr(end)/n_qpSub);
I_gauss2 = zeros(numel(n_div_arr),1);
N_arr2 = zeros(numel(n_div_arr),1);
for n_div = n_div_arr
    Xi_e_y_arr  = linspace(-1,1,n_div+1);
    Eta_e_y_arr = linspace(-1,1,n_div+1);
    xi_y = zeros(n_qpSub^2,n_div^2);
    eta_y = zeros(n_qpSub^2,n_div^2);
    
    counter = 1;
    for i_eta = 1:n_div
        Eta_e_y_sub = Eta_e_y_arr(i_eta:i_eta+1);
        for i_xi = 1:n_div
            Xi_e_y_sub = Xi_e_y_arr(i_xi:i_xi+1);
            xi_y(:,counter) = parent2ParametricSpace(Xi_e_y_sub, Q_xi);
            eta_y(:,counter) = parent2ParametricSpace(Eta_e_y_sub, Q_eta);
            counter = counter + 1;
        end
    end
    xi_y = reshape(xi_y,n_div^2*n_qpSub^2,1);
    eta_y = reshape(eta_y,n_div^2*n_qpSub^2,1);
    W2D_1 = repmat(W,n_div^2,1);
%     Xi_e_arr  = linspace(-1,1,n_div+1);
%     xi = zeros(n_qpSub,n_div);
%     counter = 1;
%     for i_xi = 1:n_div
%         Xi_e_sub = Xi_e_arr(i_xi:i_xi+1);
%         xi(:,counter) = parent2ParametricSpace(Xi_e_sub, Q1D(:,1));
%         counter = counter + 1;
%     end
%     xi = reshape(xi,n_div*n_qpSub,1);
%     W2D_1 = repmat(W1D,n_div,1);
    J_2 = 1/n_div^2;
    I_gauss2(n_div) = sum(f(xi_y,eta_y).*W2D_1*J_2);
    N_arr2(n_div) = n_div*n_qpSub;
end
    
figure(2)
totError = abs(I_gauss - I)/abs(I);
totError2 = abs(I_gauss2 - I)/abs(I);
semilogy(N_arr,totError,N_arr2,totError2)
legend('Gauss','Sub Gauss')









