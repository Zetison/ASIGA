close all
clear all %#ok

addpath ..
addpath ../utils
addpath ../models

startMatlabPool

setS135Parameters
SHBC = false;
ESBC = false;
SSBC = false;
R_i = R_o - t;
f = 30e3;
omega = 2*pi*f; % Angular frequency

defineBCstring
d_vec = [0, 0, -1];

options = struct('d_vec', d_vec,... 
                 'omega', omega, ...
                 'R_i', R_i, ...
                 'R_o', R_o, ...
                 'P_inc', P_inc, ...
                 'E', E, ...
                 'nu', nu, ...
                 'rho_s', rho_s, ...
                 'rho_f', rho_f, ...
                 'c_f', c_f);
v = -R_o(1)*d_vec;
%% Create spy matrix (requires to be in debug mode in getCoeffs.m)
% uncomment keyboard command in getCoeffs.m to get spy matrix
data = e3Dss(v, options);