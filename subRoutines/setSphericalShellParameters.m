
    
P_0 = 1; % Amplitude of incident wave
rho = 1000; % Density of fluid
rho_s = 7669; % Density of solid
rho_f_i = 1000; % Density of fluid
c = 1524;           % Speed of sound in outer (exterior) fluid domain
t = 0.008; % The thickness of the sphere
R = 5; % Midsurface of all spheres
E = 207e9; % Youngs modulus of elastic material
nu = 0.3; % Poisson ratio of elastic material
R_i = R - t/2; % Inner radius of shell
R_o = R + t/2; % Outer radius of shell
c_f_i = 340;

define_k_From_Freq = true;

if define_k_From_Freq
    f = 5e2;             % Frequency
    omega = 2*pi*f;      % Angular frequency
    k = omega/c;   % Wave number for outer fluid domain

else
    k = 2;             % Wave number for outer fluid domain
    omega = c*k;   % Angular frequency
    f = omega/(2*pi);    % Frequency
end
