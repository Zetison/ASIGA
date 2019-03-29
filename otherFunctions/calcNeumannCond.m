function dP_inc = calcNeumannCond(x_vec, n, alpha_s, beta_s, k, P_0, x_0)

k_vec = -[k*cos(beta_s)*sin(alpha_s);
          k*sin(beta_s);
          k*cos(beta_s)*cos(alpha_s)];
     
P_inc = P_0*exp(-1i*dot(k_vec,x_0))*exp(1i*dot(k_vec,x_vec));

dP_inc = 1i*dot(k_vec,n)*P_inc;