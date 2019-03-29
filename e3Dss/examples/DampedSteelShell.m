close all
clear all %#ok

addpath ..
addpath ../utils
addpath ../models

startMatlabPool

%% Define parameters
rho_f = [1000, 1000, 1.21]; % Density of fluids
rho_s = [800, 7000]; % Density of solid
c_f = [1500, 1500, 343]; % Speed of sound in fluid domains
c_s_1 = [123, 5760]; % longitudinal sound speed (compressional velocity)
c_s_2 = [33, 3130]; % transversial sound speed (shear velocity)
att = [0.1,0.01]; % Damping for youngs modulus
R_o = [0.52,0.5]; % Outer radii
epsilon = 1e-10; % solution exact when epsilon -> 0
R_i = [0.5+epsilon,0.48]; % Inner radii
d_vec = [0,0,1].'; % direction of incident wave

nFreqs = 3000; % number of frequencies in sweep

%% Calculate dependent parameters
E = c_s_2.^2.*rho_s.*(3*c_s_1.^2-4*c_s_2.^2)./(c_s_1.^2-c_s_2.^2).*(1 - 1i*att);
% E = c_s_2.^2.*rho_s.*(3*c_s_1.^2-4*c_s_2.^2)./(c_s_1.^2-c_s_2.^2);
nu = (c_s_1.^2-2*c_s_2.^2)./(2*(c_s_1.^2-c_s_2.^2));

k = linspace(10/nFreqs,10,nFreqs).'/R_o(1);
omega = k*c_f(1);   % Wave number for outer fluid domain

%% Run simulation
options = struct('d_vec', d_vec, ...
                 'omega', omega, ...
                 'R_i', R_i, ...
                 'R_o', R_o, ...
                 'E', E, ...
                 'nu', nu, ...
                 'rho_s', rho_s, ...
                 'rho_f', rho_f, ...
                 'calc_farField', true, ...
                 'c_f', c_f);
             
v = [0,0,-1];           % Compute backscattered pressure

data = e3Dss(v, options);

%% Plot results
TS = 20*log10(abs(data(1).p));
plot(k*R_o(1), TS)
set(0,'defaulttextinterpreter','latex')
xlabel('$$k_1 R_{0,1}$$')
ylabel('TS')
xlim(R_o(1)*[k(1) k(end)])



