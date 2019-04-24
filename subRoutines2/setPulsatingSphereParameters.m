rho_f = [1, NaN]; % Density of outer fluid
v_0 = -1650*1i; % The (constant) acoustic radial velocity
k = 1;
c_f = [1, NaN];
omega = c_f(1)*k;
f = omega/(2*pi);
R_o = 0.5;
t = 0.2;
P_inc = -1i*omega*rho_f(1)*v_0*R_o^2/(1-1i*k*R_o)*exp(-1i*k*R_o)*4*pi;



E = NaN;
nu = NaN;
rho_s = NaN;