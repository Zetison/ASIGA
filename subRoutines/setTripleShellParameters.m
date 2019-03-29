P_inc = 1; % Amplitude of incident wave
if 0
    rho_f = [1000 1000 1000 1000]; % Density of fluids
    c_f = [1500, 1500, 1500, 1500]; % Speed of sound in fluid domains
else
    rho_f = [1000 1000 1000 0.01]; % Density of fluids
    c_f = [1500 1500 1500 0.3]; % Speed of sound in fluid domains
end
rho_s = 7850*[1,1,1]; % Density of solid
E = 210e9*[1,1,1]; % Youngs modulus of elastic material
% E = 210e6*[1,1,1]; % Youngs modulus of elastic material
nu = 0.3*[1,1,1]; % Poisson ratio of elastic material
t = [0.2, 0.8, 0.1];
R_o = [6, 4, 2];