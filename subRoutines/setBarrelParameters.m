P_inc = 1; % Amplitude of incident wave
rho_s = 7850; % Density of solid
rho_f = [1000, 1000]; % Density of fluids
c_f = [1500, 1.2];  % Speed of sound in fluid domains

t = 0.2; % The thickness of the sphere
E = 210e9; % Youngs modulus of elastic material
nu = 0.3; % Poisson ratio of elastic material
R_o = 2.5; % Outer radius of shell
R_i = R_o - t; % Inner radius of shell
L = 20;

f = 1e2;             % Frequency
omega = 2*pi*f;      % Angular frequency
k = omega/c;   % Wave number for outer fluid domain

alpha_s = 90*pi/180;
beta_s = 0*pi/180;    