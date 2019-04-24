
P_0 = 1; % Amplitude of incident wave
rho_f = 1000; % Density of fluid
rho_s = 7669; % Density of solid
c_f = 1524; % Speed of sound in fluid
R = 5; % Midsurface of all spheres
t = 0.15; % The thickness of the sphere
E = 207e9; % Youngs modulus of elastic material
nu = 0.3; % Poisson ratio of elastic material
R_i = R - t/2; % Inner radius of shell
R_o = R + t/2; % Outer radius of shell

% Define Lam� parameters
lambda = nu*E/((1+nu)*(1-2*nu));
mu = E/(2*(1+nu));

c_s = sqrt(E/((1-nu^2)*rho_s)); % Speed of sound in elastic material ??????????????????????????????????????
c_1 = sqrt((lambda+2*mu)/rho_s);
c_2 = sqrt(mu/rho_s);