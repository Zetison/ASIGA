% This script produces Figure 19 in 

close all
clear all %#ok

addpath ..
addpath ../utils
addpath ../models

% pathToResults = '../../../results/e3Dss/';
pathToResults = '../results';

startMatlabPool

model = 'S15';
SHBC = 0;
SSBC = 0;
ESBC = 0; 
f_c = 1500;
% applyLoad = 'pointCharge';
applyLoad = 'radialPulsation';
switch applyLoad
    case 'pointCharge'
        ESBC = 1;
        d_vec = -[-sqrt(r_s^2-R_a^2), R_a, 0]';  
        d_vec = d_vec/norm(d_vec); 
        T = 120/f_c;
        N = 2^10;
    case 'planeWave'
        d_vec = [1, 0, 0]';   
        ESBC = 1;
        T = 120/f_c;
        N = 2^10;
    case 'radialPulsation'
        d_vec = [1, 0, 0]';  
        SHBC = 1; 
        T = 120/f_c;
        N = 2^11;
end
setS15Parameters

R_i = R_o - t;
defineBCstring
R_a = 1.5*R_o(1);
r_s = 2*R_o(1);

options = struct('d_vec', d_vec, ...
                 'R_i', R_i, ...
                 'R_o', R_o, ...
                 'P_inc', P_inc, ...
                 'E', E, ...
                 'nu', nu, ...
                 'rho_s', rho_s, ...
                 'rho_f', rho_f, ...
                 'c_f', c_f, ...
                 'applyLoad', applyLoad, ...
                 'plotTimeOscillation', 0, ...
                 'plotInTimeDomain', 1, ...
                 'SHBC', SHBC, ...
                 'ESBC', ESBC, ...
                 'SSBC', SSBC, ...
                 'R_a', R_a, ...
                 'r_s', r_s, ...
                 'f_c', f_c, ...
                 'N', N, ...
                 'T', T,...
                 'computeForSolidDomain', 0);

options.Eps = 1e-8;
extraPts = 30;
folderName = [pathToResults '/paraviewResults/' model '_' BC '/'];
if ~exist(folderName, 'dir')
    mkdir(folderName);
end

createParaviewFiles_e3Dss(extraPts, folderName, options)     