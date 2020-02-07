P_inc = 1; % Amplitude of incident wave
rho_s = 7850; % Density of solid
rho_f = [1000, 1000]; % Density of fluids
c_f = [1500, 1500];  % Speed of sound in fluid domains
E = 210e9; % Youngs modulus of elastic material
nu = 0.3; % Poisson ratio of elastic material


t = 0.01;   % Thickness

a = 7.0;
b = 3.5;
L = 42.0;

g2 = 6.5;
g3 = 6.5;
alpha = 18*pi/180;
beta = 240*pi/180;

c = 4.0;
s = 1.2;

b_ls = 2*s;
b_us = 2;
l_ls = 13;
l_us = 12.3;
h_s = 3.5;
delta_s = 0.2;
x_s = a-19;

l_lm = 2.6;
b_lm = 0.4;
l_um = 2.35;
b_um = 0.3;
h_m = 3.5;
delta_m = 0.25;
x_m = -51.9;


C_4 = c + b*cos(beta/2);
C_3 = b*sin(beta/2)-s;
C_2 = (2*C_4+C_3*tan(beta/2))/C_3^3;
C_1 = -(3*C_4+C_3*tan(beta/2))/C_3^2;

nPanels = 5;
dyPanels = C_3/nPanels;
P_panels = @(y) c + C_1*(y-s).^2 + C_2*(y-s).^3;

l_ld = 2.6;
b_ld = 2*(c-P_panels(s+dyPanels));
l_ud = 2.35;
b_ud = 0.22;
h_d = b-s;
delta_d = 0.25;
x_d = -4;
