P_0 = 1; % Amplitude of incident wave
rho = 1000; % Density of fluid
rho_f_i = 1.2; % Density of fluid inner part
rho_s = 7850; % Density of solid
c = 1524; % Speed of sound in outer fluid
c_f_i = 340;  % Speed of sound in inner fluid

t = 0.5; % The thickness of the sphere
E = 210e9; % Youngs modulus of elastic material
nu = 0.3; % Poisson ratio of elastic material
R_o = 3; % Outer radius of shell
R_i = R_o - t; % Inner radius of shell
L = 13; %40
% L = 6;

if true
    f = 5e2;             % Frequency
    omega = 2*pi*f;      % Angular frequency
    k_wn_o = omega/c;   % Wave number for outer fluid domain

else
    k_wn_o = 1;             % Wave number for outer fluid domain
    omega = c*k_wn_o;   % Angular frequency
    f = omega/(2*pi);    % Frequency
end
