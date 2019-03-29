P_inc = 1; % Amplitude of incident wave
rho_s = 7669; % Density of solid
rho_f = [1000, 1000]; % Density of fluids
c_f = [1524, 1524];  % Speed of sound in fluid domains
t = 0.15; % The thickness of the sphere
E = 207e9; % Youngs modulus of elastic material
nu = 0.3; % Poisson ratio of elastic material

R = 5; % Midsurface radius
R_o = R+t/2; % Outer radius of shell