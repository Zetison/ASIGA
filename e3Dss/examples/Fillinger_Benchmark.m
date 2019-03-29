close all
clear all %#ok

addpath ..
addpath ../utils
addpath ../models

alpha_s = 0;
beta_s = 0;

beta_f = beta_s;
alpha_f = alpha_s;

d_vec = -[cos(beta_s)*cos(alpha_s);
          cos(beta_s)*sin(alpha_s);
          sin(beta_s)];

P_inc = 1; % Amplitude of incident wave
rho_f = 1000;
c_f = 1500;
k = 10.^linspace(-1,3,3000);
R_o = 1; % Outer radius of shell

omega = c_f*k; % Angular frequency

options = struct('d_vec', d_vec,... 
                 'omega', omega, ...
                 'P_inc', P_inc, ...
                 'rho_f', rho_f, ...
                 'calc_farField', 1, ...
                 'c_f', c_f);
             
             
             
v = R_o(1)*[cos(beta_f)*cos(alpha_f); cos(beta_f)*sin(alpha_f); sin(beta_f)*ones(size(alpha_f))]';
data = e3Dss(v, options);
TS_inf = 20*log10(R_o/2)*ones(size(k));
TS = 20*log10(abs(data(1).p));

% figure(1)
% semilogx(k, TS,'DisplayName','Analytic')
% figure(2)
% f = omega/(2*pi);
% error_pAbs = 100*abs(abs(data(1).p)-abs(R_o/2))./abs(data(1).p);
% loglog(f,error_pAbs)
% figure(3)
% error_p = 100*abs(data(1).p-P_inc*R_o/2*exp(-2*1i*k*R_o))./abs(data(1).p);
% loglog(f,error_p)
% % printResultsToFile2('../../results/_studies/Fillinger/asymptotic', f.', TS_inf.', 'f','TS');
% % printResultsToFile2('../../results/_studies/Fillinger/asymptotic', f.', error_pAbs.', 'f','error_pAbs');
% % printResultsToFile2('../../results/_studies/Fillinger/asymptotic', f.', error_p.', 'f','error_p');
% return
figure(1)
semilogx(k, TS-TS_inf,'DisplayName','Analytic')
ylim([-30,10])
hold on
k = 10.^linspace(-1,-0.35,2);
TS = 20*log10(5/6*k.^2*R_o^3);
semilogx(k, TS-TS_inf,'green','DisplayName','Asymptotic')
legend('off');
legend('show');

k = 10.^linspace(0.2,2,2);
TS = 10*log10(R_o^2/4)*ones(size(k));
semilogx(k, TS-TS_inf,'green','DisplayName','Asymptotic')


xlabel('$$k\cdot a$$','interpreter','latex')
ylabel('$$\mathrm{TS}-\mathrm{TS}_{\infty}$$','interpreter','latex')
