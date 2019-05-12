close all
% 
% rule = 'Freud';
type = 'mp';
% type = 'double';
switch type
    case 'single'
        Eps = 1e-7;
        PI = pi;
        O = 0;
    case 'double'
        Eps = eps;
        PI = pi;
        O = 0;
    case 'sym'
        Eps = 1e-100;
        digits(2000);
        PI = vpa('pi');
        O = vpa('0');
    case 'mp'
        Eps = 1e-100;
        mp.Digits(100);
        PI = mp('pi');
        O = mp('0');
end
rule = 'Legendre';
N = 143;
I_gauss = zeros(N,1,type);
N_arr = 1:N;
omega = mp('100');
f = @(x) cos(omega*x);
a = 0.3;
% f = @(x) abs(x-a);
c = 5;
% f = @(x) abs(x-a).*(x-a).^c;
% f = @(x) exp(x).*sin(x);
% f = @(x) exp(1i*x+x).*legendre_(7,x);
for n = N_arr
    n
    switch rule
        case 'Legendre'
            [Q1D,W1D] = getQuadFromFile(n,type);
%             [Q1D,W1D]=lgwt(n);
        case 'Freud'
            r = 4;
            [Q1D,W1D] = gaussFreudQuad(n,r); 
    end

    I_gauss(n) = sum(f(Q1D).*W1D);
end
% switch rule
%     case 'Legendre'
%         I = integral(f,-1,1);
%     case 'Freud'
%         integrand = @(x) f(x).*exp(-x.^(r+1));
%         I = integral(integrand,0,inf);
% end
I = 2*sin(omega)/omega;
figure(1)
totError = abs(I_gauss - I)
% semilogy(N_arr,totError)
% loglog(N_arr,totError,N_arr,N_arr.^(-(c+2)))
% semilogy(N_arr,factorial(N_arr).^2./factorial(2*N_arr))
% hold on 
semilogy(N_arr,totError,'DisplayName','Gaussian quadrature error: $$\left|\int_{-1}^1 f(x)\,\mathrm{d}x - \sum_{i=1}^n f(x_i)w_i\right|$$')
hold on
N_arr = mp(N_arr);
semilogy(N_arr,2.^(2*N_arr+1).*factorial(N_arr).^4./(2*N_arr+1)./factorial(2*N_arr).^3.*omega.^(2*N_arr),'DisplayName', ...
    'Upper bound on error: $$\frac{2^{2n+1}(n!)^4}{(2n+1)[(2n)!]^3}\max_{\xi\in[-1,1]}f^{(2n)}(\xi)$$')
% legend('error','$$2.^(2*N_arr+1).*factorial(N_arr).^4./(2*N_arr+1)./factorial(2*N_arr).^3.*(pi*k).^(2*N_arr)$$')
leg1 = legend('show','Location','northwest');
set(leg1,'Interpreter','latex');
xlabel('$n$','interpreter','latex')
% disp(rms(totError))
title(['Gaussian integration of $f(x) = \cos(\omega x),\, \omega =$ ', num2str(double(omega))],'interpreter','latex')
hold off
savefig('../gaussErrorHP')
return
% for type = {'double','mp'}
%     noqpMax = 1000;
%     W = cell(noqpMax,1);
%     Q = cell(noqpMax,1);
%     mp.Digits(100);
%     for ii = 1:noqpMax
%         [Q{ii},W{ii}] = getQuadFromFile(ii,type{1});
%     end
%     quadData.Q = Q;
%     quadData.W = W;
%     save(['integration/quadData_' type{1}],'quadData')
% end

% close all
% 
% for errorExp = 10
%     rule = 'Legendre';
%     N = 64;
%     n_arr = 1:N;
%     pik_gauss = zeros(size(n_arr));
%     k = linspace(0,100,1000);
%     f = @(x) cos(pi*x*k);
%     % f = @(x) cos(100*x);
%     % f = @(x) exp(x).*sin(x);
%     for i = 1:numel(n_arr)
%         n = n_arr(i);
%         [W1D,Q1D] = gaussianQuadNURBS(n); 
%         I = 2*sin(pi*k)./k/pi;
%         I_gauss = sum(W1D.'*f(Q1D),1);
%         totError = abs(I_gauss - I)./abs(I);
%         temp = k(min(find(totError > 10^(-errorExp))));
%         pik_gauss(i) = temp*pi;
% 
%     end
%     figure(1)
%     plot(pik_gauss, n_arr,'DisplayName',num2str(errorExp))
%     xlabel k
%     ylabel n
%     hold on 
% end
% legend('off');
% legend('show','Location','best');

% figure(2)
% pik = linspace(0,100,1000);
% totError = zeros(size(pik));
% for i = 1:numel(pik)
%     k = pik(i);
%     f = @(x) cos(x*k);
%     n = ceil(-7.7./(1+k)+0.7*k+10);
%     [W1D,Q1D] = gaussianQuadNURBS(n); 
%     I = 2*sin(k)/k;
%     I_gauss = W1D.'*f(Q1D);
%     totError(i) = abs(I_gauss - I)/abs(I);
% end
% semilogy(pik,totError)
%     
% 
% 
