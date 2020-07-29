function container = setBeTSSi_M31Parameters()

P_inc = 1; % Amplitude of incident wave
rho_s = 7850; % Density of solid
rho_f = [1000, 1000, 1.2]; % Density of fluids
c_f = [1500, 1500, 340];  % Speed of sound in fluid domains
E = 210e9; % Youngs modulus of elastic material
nu = 0.3; % Poisson ratio of elastic material


t = [0.008,0.02];   % Thickness
R_o1 = 5; % Outer radius
R_o2 = 3; % Outer radius
R_i1 = R_o1-t(1); % Inner radius
R_i2 = R_o2-t(1); % Inner radius

R_o = 3;
R_i = R_o-t(2);
L = [49-R_o1-R_o2,43-R_o];

putVariablesIntoContainer