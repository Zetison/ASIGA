function container = setTAPParameters()

% elastic parameters from nickel metal
P_inc = 1; % Amplitude of incident wave
rho_s = 8800; % Density of solid
rho_f = [1000, 1.2]; % Density of fluids
c_f = [1510, 340];  % Speed of sound in fluid domains
E = 210e9; % Youngs modulus of elastic material
nu = 0.3; % Poisson ratio of elastic material

t = 0.023; % BeTSSi thickness
% L = 120/1.0936; % Shirron example 2
L = 42.3; % Shirron example 2

% R_o = 5.5/1.0936;
R_o = 3.09;

putVariablesIntoContainer