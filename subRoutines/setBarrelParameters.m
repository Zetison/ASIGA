P_inc = 1; % Amplitude of incident wave
rho_s = 7850; % Density of solid
rho_f = [1000, 340]; % Density of fluids
c_f = [1500, 1.2];  % Speed of sound in fluid domains

t = 0.1; % The thickness of the barrel
E = 210e9; % Youngs modulus of elastic material
nu = 0.3; % Poisson ratio of elastic material
R_o = 1; % Outer radius of shell
R_i = R_o - t; % Inner radius of shell
L = 4*pi/4;