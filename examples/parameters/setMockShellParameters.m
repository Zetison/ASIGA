function container = setMockShellParameters()

P_inc = 1; % Amplitude of incident wave
rho_s = 7850; % Density of solid
rho_f = [1000, 1000]; % Density of fluids
c_f = [1500, 1500];  % Speed of sound in fluid domains
E = 210e9; % Youngs modulus of elastic material
nu = 0.3; % Poisson ratio of elastic material

R_o = 1;
% t = 0.0406; % Thickness of the Mock shell of Ihlenburg
t = 0.008; % BeTSSi thickness
% L = 4*R_o; % Shirron example 1
L = 10*R_o; % Shirron example 2

% mult = 1;

putVariablesIntoContainer