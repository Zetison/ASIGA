close all
clear all %#ok

addpath ..
addpath ../utils
addpath ../models

pathToResults = '../../plotData/e3Dss/';
% pathToResults = pathToResults '';

startMatlabPool

%% Create plot of bessel functinos for large n
xi = 500;
type = 2;
switch type
    case 1
        N = 550;
    case 2
        N = 1050;
end
npts = N;
B = zeros(npts,1);
B1 = zeros(npts,1);
B2 = zeros(npts,1);
N_arr = linspace(1,N,npts).';
parfor i = 1:length(N_arr)
    n = N_arr(i);
    B1(i) = bessel_s(n,xi,1);
    B2(i) = bessel_s(n,xi,2);
end
semilogy(N_arr,abs(B1),N_arr,abs(B2))
xlim([0,N])
xlabel('$n$','interpreter','latex')
ylabel('Magnitude of Bessel functions')
switch type
    case 1
        ylim([2e-10,5e3])
        savefig([pathToResults 'Figure4a'])
    case 2
        ylim([1e-223,1e223])
        savefig([pathToResults 'Figure4b'])
end
legend({'$$|\mathrm{j}_n(500)|$$','$$|\mathrm{y}_n(500)|$$'},'interpreter','latex','location','northwest')

% switch type
%     case 1
%         ylim([1e-10 1e4])
%         printResultsToFile(pathToResults 'besseljForLargeN1', N_arr, abs(B1), [], 0, 1)
%         printResultsToFile(pathToResults 'besselyForLargeN1', N_arr, abs(B2), [], 0, 1)
%     case 2
%         ylim([1e-220 1e220])
%         printResultsToFile(pathToResults 'besseljForLargeN2', N_arr, abs(B1), [], 0, 1)
%         printResultsToFile(pathToResults 'besselyForLargeN2', N_arr, abs(B2), [], 0, 1)
% end