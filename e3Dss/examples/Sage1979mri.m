close all
clear all %#ok

addpath ..
addpath ../utils
addpath ../models

startMatlabPool

%% Define parameters
rho_f = [1025, 1.293]; % Density of fluids
rho_s = NaN; % Density of solid
c_f = [1531, 346.2]; % Speed of sound in fluid domains
R_o = 1; % Outer radii
R_i = []; % Inner radii
E = NaN;
nu = NaN;

E = 2e6;
nu = 0.46;
rho_s = 1; % Density of solid
R_i = R_o-1e-6; % Inner radii
% rho_f = rho_f(1);
% c_f = c_f(1);
% rho_s = [];
nFreqs = 3000;

%% Calculate dependent parameters

k = [10.^linspace(-3,0,nFreqs)'; linspace(1+5/nFreqs,5,nFreqs)'];
omega = k*c_f(1);

d_vec = [0,0,1].';

%%%%%%%%%
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
             
omega = k*c_f(1);   % Wave number for outer fluid domain
options.omega = omega;
v = [0,0,-1];           % Compute backscattered pressure

data = e3Dss(v, options);

%% Plot results
SPL_Sage1 = importdata('../models/Sage1979mri/Figure1.csv');
SPL_Sage2 = importdata('../models/Sage1979mri/Figure4.csv');
SPL_Sage = [SPL_Sage1(SPL_Sage1(:,1) < 1,1), SPL_Sage1(SPL_Sage1(:,1) < 1,2);
            SPL_Sage2(SPL_Sage2(:,1) > 1,1), SPL_Sage2(SPL_Sage2(:,1) > 1,2)];
        
% semilogy(k*R_o, 4*abs(data(1).p).^2/R_o^2/pi,SPL_Sage(:,1),SPL_Sage(:,2))
loglog(k*R_o, 4*abs(data(1).p).^2/R_o^2/pi,SPL_Sage(:,1),SPL_Sage(:,2))
warning('y axis scaled with another pi which is not consistent with Sage1979mri...')
set(0,'defaulttextinterpreter','latex')
legend({'Present work', 'Reference Solution from Sage (1979)'})
xlabel('$$k_1 R_{0,1}$$')
ylabel('$$\frac{\sigma}{\pi R_{o,1}^2}$$')
xlim(R_o*[k(1) k(end)])



