function container = setBeTSSi_M3Parameters()

P_inc = 1; % Amplitude of incident wave
rho_s = 7850; % Density of solid
rho_f = [1000, 1000]; % Density of fluids
c_f = [1500, 1500];  % Speed of sound in fluid domains
E = 210e9; % Youngs modulus of elastic material
nu = 0.3; % Poisson ratio of elastic material


t = 0.008;   % Thickness
R_o1 = 5; % Outer radius
R_o2 = 3; % Outer radius
R_i1 = R_o1-t; % Inner radius
R_i2 = R_o2-t; % Inner radius
% R_o1 = 5 + t/2; % Outer radius
% R_o2 = 3 + t/2; % Outer radius
% R_i1 = 5-t/2; % Inner radius
% R_i2 = 3-t/2; % Inner radius
L = 49-R_o1-R_o2;

putVariablesIntoContainer