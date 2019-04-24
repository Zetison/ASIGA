k_wn_o = 3;
c_f = 1500;
omega = k_wn_o*c_f;
k_vec = k_wn_o*[1; 2]/sqrt(5);

u_analytic = @(x) P_0*exp(1i*dot2(x,k_vec));
Grad = @(x) 1i*k_vec*u_analytic(x);