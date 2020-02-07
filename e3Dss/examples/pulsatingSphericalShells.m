close all
clear all %#ok
% 
addpath ..
addpath ../utils
addpath ../models

%% Test benchmark models

beta_f = pi/2;
alpha_f = 0;
npts = 1000;

f = 1e3; % frequencies in Hertz
model = 'S15';
% model = 'S1';
% applyLoad = 'planeWave';
applyLoad = 'radialPulsation';

switch model
    case 'S1'
        setS1Parameters
        R_i = R_o - t;
        r_arr1 = linspace(R_o(1),2*R_o(1),npts).';
        r_arr2 = linspace(R_i(1),R_o(1),npts).';
        r_arr3 = linspace(0,R_i(1),npts).';
    case 'S15'
        setS15Parameters
        R_i = R_o - t;
        r_arr1 = linspace(R_o(1),2*R_o(1),npts).';
        r_arr2 = linspace(R_i(1),R_o(1),npts).';
        r_arr3 = linspace(R_o(2),R_i(1),npts).';
end
v1 = r_arr1*[cos(beta_f)*cos(alpha_f), cos(beta_f)*sin(alpha_f), sin(beta_f)*ones(size(alpha_f))];
v2 = r_arr2*[cos(beta_f)*cos(alpha_f), cos(beta_f)*sin(alpha_f), sin(beta_f)*ones(size(alpha_f))];
v3 = r_arr3*[cos(beta_f)*cos(alpha_f), cos(beta_f)*sin(alpha_f), sin(beta_f)*ones(size(alpha_f))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SHBC
switch model
    case 'S1'
        setS1Parameters
    case 'S15'
        setS15Parameters
end
SHBC = true;
ESBC = false;
SSBC = false;
omega = 2*pi*f; % Angular frequency

defineBCstring

k = omega./c_f;
options = struct('omega', omega, ...
                 'R_i', R_i, ...
                 'R_o', R_o, ...
                 'P_inc', P_inc, ...
                 'E', E, ...
                 'nu', nu, ...
                 'rho_s', rho_s, ...
                 'rho_f', rho_f, ...
                 'calc_farField', 0, ...
                 'calc_sigma_rr', 1, ...
                 'applyLoad',applyLoad,...
                 'c_f', c_f);
switch applyLoad
    case 'planeWave'
        d_vec = [0;0;1];
        p_inc = @(r) P_inc*exp(1i*k(1)*r);
    case 'radialPulsation'
        p_inc = @(r) P_inc*R_o(1)*exp(-1i*k(1)*(r-R_o(1)))./r;
end
% data = e3Dss(v1, options);
data = e3Dss({v1, v2, v3}, options);
figure(42)
plot([r_arr3; r_arr2; r_arr1], [real(data(2).p); -real(data(1).sigma_rr); real(p_inc(r_arr1) + data(1).p)])
% plot([r_arr2; r_arr1], [-real(data(1).sigma_rr); real(p_inc(r_arr1) + data(1).p)])
% plot(r_arr1, real(p_inc(r_arr1) + data(1).p))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ESBC
switch model
    case 'S1'
        setS1Parameters
    case 'S15'
        setS15Parameters
end
SHBC = false;
ESBC = false;
SSBC = true;
R_i = R_o - t;
omega = 2*pi*f; % Angular frequency

defineBCstring

k = omega./c_f;
options = struct('omega', omega, ...
                 'R_i', R_i, ...
                 'R_o', R_o, ...
                 'P_inc', P_inc, ...
                 'E', E, ...
                 'nu', nu, ...
                 'rho_s', rho_s, ...
                 'rho_f', rho_f, ...
                 'calc_farField', 0, ...
                 'calc_sigma_rr', 1, ...
                 'applyLoad',applyLoad,...
                 'c_f', c_f);
% data = e3Dss({v1, v2}, options);
data = e3Dss({v1, v2, v3}, options);
hold on
% plot(r_arr1, real(data(1).p),r_arr2, real(data(2).p))
% plot([r_arr2; r_arr1], [-real(data(1).sigma_rr); real(p_inc(r_arr1) + data(1).p)])
plot([r_arr3; r_arr2; r_arr1], [real(data(2).p); -real(data(1).sigma_rr); real(p_inc(r_arr1) + data(1).p)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NNBC
switch model
    case 'S1'
        setS1Parameters
    case 'S15'
        setS15Parameters
end
SHBC = false;
ESBC = false;
SSBC = false;
R_i = R_o - t;
omega = 2*pi*f; % Angular frequency

defineBCstring

k = omega./c_f;
options = struct('omega', omega, ...
                 'R_i', R_i, ...
                 'R_o', R_o, ...
                 'P_inc', P_inc, ...
                 'E', E, ...
                 'nu', nu, ...
                 'rho_s', rho_s, ...
                 'rho_f', rho_f, ...
                 'calc_farField', 0, ...
                 'calc_sigma_rr', 1, ...
                 'applyLoad',applyLoad,...
                 'c_f', c_f);
data = e3Dss({v1, v2, v3}, options);
hold on
% plot(r_arr1, real(data(1).p),r_arr2, real(data(2).p))
plot([r_arr3; r_arr2; r_arr1], [real(data(2).p); -real(data(1).sigma_rr); real(p_inc(r_arr1) + data(1).p)])
legend('SHBC','SSBC','NNBC')
legend show
hold off

