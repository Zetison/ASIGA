close all
clear all %#ok

addpath ..
addpath ../utils
addpath ../models

pathToResults = '../../../results/e3Dss/';
% pathToResults = '../results';

P_inc = 1; % Amplitude of incident wave
rho_f = 1000; % Density of outer fluid
rho_s = 7850; % Density of solid
c_f = 1500;   % Speed of sound in outer (exterior) fluid domain
t = 0.05; % The thickness of the sphere
R_o = 1; % Outer radius of shell
R_i = R_o-t; % Inner radius of shell
E = 2.0e11; % Youngs modulus of elastic material
nu = 0.3; % Poisson ratio of elastic material

if false
    f = (1430:12:4290)';
    omega = 2*pi*f;   % Wave number for outer fluid domain
    k = omega/c_f(1);
else
    k = (6:0.05:18)';
    omega = k*c_f(1);
    f = omega/(2*pi);
end
d_vec = [0,0,1].';
options = struct('d_vec', d_vec, ...
                 'omega', omega, ...
                 'R_i', R_i, ...
                 'R_o', R_o, ...
                 'P_inc', P_inc, ...
                 'E', E, ...
                 'nu', nu, ...
                 'rho_s', rho_s, ...
                 'rho_f', rho_f, ...
                 'c_f', c_f);
             
options.omega = omega;

v = R_o*[0,0,-1];
data = e3Dss(v, options);

figure(3)
real_p_Hetmaniuk = importdata('../models/Hetmaniuk2012raa/Figure17.csv');
p_inc = @(v) exp(1i*dot3(v,d_vec)*k.');
plot(f, real(data(1).p),'DisplayName','Exact')
title('Figure 17 in Hetmaniuk2012raa')
hold on
plot(real_p_Hetmaniuk(:,1), real_p_Hetmaniuk(:,2),'DisplayName','Hetmaniuk2012raa')
xlabel('Frequency [Hz]')
xlim([f(1), f(end)])
ylim([-15, 15])
ylabel('Real part of pressure')  
legend('show');



