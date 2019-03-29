close all
% 
% rule = 'Freud';
rule = 'Legendre';
N = 63;
I_gauss = zeros(N,1);
I_gauss2 = zeros(N,1);
N_arr = 1:N;
N_arr2 = zeros(size(N_arr));
epsilon = 0.1;
f = @(x,y) 1./sqrt(epsilon^2+x.^2+y.^2);
NN = 5;
[W2D2,Q2D2] = gaussianQuadNURBS(NN,NN);
% f = @(x) exp(x).*sin(x);
% f = @(x) exp(1i*x+x).*legendre_(7,x);
for n = N_arr
    [W2D,Q2D] = gaussianQuadNURBS(n,n);

    I_gauss(n) = sum(f(Q2D(:,1),Q2D(:,2)).*W2D);
    
    Xi_e_arr  = linspace(-1,1,n+1);
    Eta_e_arr = linspace(-1,1,n+1);
    for i_eta = 1:n
        Eta_e_sub = Eta_e_arr(i_eta:i_eta+1);
        for i_xi = 1:n
            Xi_e_sub = Xi_e_arr(i_xi:i_xi+1);

            x  = parent2ParametricSpace(Xi_e_sub, Q2D2(:,1));
            y = parent2ParametricSpace(Eta_e_sub,Q2D2(:,2));
            
            J_2 = 0.25*(Xi_e_sub(2)-Xi_e_sub(1))*(Eta_e_sub(2)-Eta_e_sub(1));
            I_gauss2(n) = I_gauss2(n) + sum(f(x,y).*W2D2.*J_2);
        end
    end
    N_arr2(n) = NN*n;
end
I = integral2(f,-1,1,-1,1,'AbsTol',1e-14,'RelTol',1e-14);
% I = 2*sin(k)/k;
figure(1)
totError = abs(I_gauss - I)/abs(I)
totError2 = abs(I_gauss2 - I)/abs(I)
loglog(N_arr.^2,totError,N_arr2.^2,totError2)
% loglog(N_arr,totError,N_arr,N_arr.^(-(c+2)))
hold on 
% semilogy(N_arr,factorial(N_arr).^2./factorial(2*N_arr))
% semilogy(N_arr,2.^(2*N_arr+1).*factorial(N_arr).^4./(2*N_arr+1)./factorial(2*N_arr).^3.*k.^(2*N_arr))
legend('error','error2')
% disp(rms(totError))
hold off
