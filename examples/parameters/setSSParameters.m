function layer = setSSParameters(noDomains)

rho_s = 7850; % Density of solid
rho_f = 1000; % Density of fluids
c_f = 1500;  % Speed of sound in fluid domains
t = 0.008; % The thickness of the sphere
E = 210e9; % Youngs modulus of elastic material
nu = 0.3; % Poisson ratio of elastic material

R_o = 0.5;

convertToLayerFormat
layer = layer(1:noDomains);