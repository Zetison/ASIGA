close all
clear all %#ok

addpath ..
addpath ../utils
addpath ../models

c_f = 1500;
k = 20;
omega = k*c_f; % Angular frequency
R_o = 1;


options = struct('d_vec', [0;0;1].', ... 
                 'omega', omega, ...
                 'R_o', R_o, ...
                 'c_f', c_f);
             
specialValues = [];
theta_arr = linspace(0,pi/2,1000);
% theta_arr = linspace(0,2*pi,1000);
v = 22/k*[sin(theta_arr); 0*zeros(size(theta_arr)); cos(theta_arr)]';
data = e3Dss(v, options);

p_inc = @(v) exp(1i*k*v(:,3));
plot(theta_arr, 20*log10(abs(data(1).p+p_inc(v))))

legend('Abs total field dB')
xlabel('$$\theta$$','interpreter','latex')
