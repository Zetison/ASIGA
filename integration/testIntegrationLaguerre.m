clear all
close all
r = 5;
alpha = -r/(1+r);

integrand = @(q) cos(q);
integrand = @(q) 1/(r+1)*q.^(5/(r+1));
integrand = @(q) 1/(r+1)*q.^(5/(r+1));

N = 30;
N_arr = 1:N;
I_gauss = zeros(N,1);
for n = N_arr
    [Q1D, W1D] = gaussLaguerreQuad(n, r);
    gaussApprox = 0;

    for i = 1:n
        w = W1D(i);
        q = Q1D(i);
        gaussApprox = gaussApprox + w*integrand(q);
    end
    I_gauss(n) = gaussApprox;
end
I = [0.5, 1.376996331853153438664376624047155467474, 2.305339340080741214795588515899745485588, ...
     3.260815786952047287241917042250172207226, 4.230672773988538325039644285926158316913, 5.208954824504981875328948024075465843530];
I = [120, 1, 1/3, .2215567313631895034122709354176431478498, .1836337484799521221281903310371660801373, 1/6];
totError = abs(I_gauss - I(r+1))./abs(I(r+1));
semilogy(N_arr,totError,'*-')
hold on

semilogy(N_arr,factorial(N_arr).*gamma(N_arr+alpha+1)./factorial(2*N_arr))
% disp(rms(totError))