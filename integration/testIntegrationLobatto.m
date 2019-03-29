close all
% 
rule = 'Legendre';
% rule = 'Lobatto';
N = 20;
I_gauss = zeros(N,1);
N_arr = 1:N;
p = N-3;
f = @(x) x.^p + x.^(p-1);
for n = N_arr
    n
    switch rule
        case 'Legendre'
            [W1D,Q1D] = gaussianQuadNURBS(n); 
        case 'Lobatto'
            [Q1D, W1D] = GLLpoints(n);
    end
    I_gauss(n) = sum(f(Q1D).*W1D);
end
I = 1/(p+1)+1/p-(-1)^(p+1)/(p+1)-(-1)^p/p;
figure(1)
totError = abs(I_gauss - I)/abs(I)
semilogy(N_arr,totError)
hold on 
% semilogy(N_arr,factorial(N_arr).^2./factorial(2*N_arr))
% semilogy(N_arr,2.^(2*N_arr+1).*factorial(N_arr).^4./(2*N_arr+1)./factorial(2*N_arr).^3.*k.^(2*N_arr))
legend('error','refline')
% disp(rms(totError))
hold off
