
% u_analytic = @(r,theta,phi) 1/E*p_i*R_i^2/(R_o^2-R_i^2)*((1+nu)*(1-2*nu)*r + R_o^2*(1+nu)/r);

coeffs = [(lambda+2*mu)*sqrt(rho_s)*omega*bessel1Deriv(1,sqrt(rho_s)*omega*R_o) + 2*lambda/R_o*bessel1(1,sqrt(rho_s)*omega*R_o), ...
                (lambda+2*mu)*sqrt(rho_s)*omega*bessel2Deriv(1,sqrt(rho_s)*omega*R_o) + 2*lambda/R_o*bessel2(1,sqrt(rho_s)*omega*R_o);
          (lambda+2*mu)*sqrt(rho_s)*omega*bessel1Deriv(1,sqrt(rho_s)*omega*R_i) + 2*lambda/R_i*bessel1(1,sqrt(rho_s)*omega*R_i), ...
                (lambda+2*mu)*sqrt(rho_s)*omega*bessel2Deriv(1,sqrt(rho_s)*omega*R_i) + 2*lambda/R_i*bessel2(1,sqrt(rho_s)*omega*R_i)]\[-p_o; -p_i];
C_1 = coeffs(1);
C_2 = coeffs(2);

u_analytic = @(r,theta,phi) C_1*bessel1(1,sqrt(rho_s)*omega*r) + C_2*bessel2(1,sqrt(rho_s)*omega*r);

sigma_rr = @(r,theta,phi) (lambda+2*mu)*sqrt(rho_s)*omega*(C_1*bessel1Deriv(1,sqrt(rho_s)*omega*r) + C_2*bessel2Deriv(1,sqrt(rho_s)*omega*r)) ...
                                + 2*lambda*(C_1*bessel1(1,sqrt(rho_s)*omega*r) + C_2*bessel2(1,sqrt(rho_s)*omega*r))/r;
sigma_thetatheta = @(r,theta,phi) lambda*sqrt(rho_s)*omega*(C_1*bessel1Deriv(1,sqrt(rho_s)*omega*r) + C_2*bessel2Deriv(1,sqrt(rho_s)*omega*r)) ...
                                    + (2*lambda+2*mu)*(C_1*bessel1(1,sqrt(rho_s)*omega*r) + C_2*bessel2(1,sqrt(rho_s)*omega*r))/r;
sigma_phiphi = @(r,theta,phi) lambda*sqrt(rho_s)*omega*(C_1*bessel1Deriv(1,sqrt(rho_s)*omega*r) + C_2*bessel2Deriv(1,sqrt(rho_s)*omega*r)) ...
                                    + (2*lambda+2*mu)*(C_1*bessel1(1,sqrt(rho_s)*omega*r) + C_2*bessel2(1,sqrt(rho_s)*omega*r))/r;
sigma_thetaphi = @(r,theta,phi) 0;
sigma_rphi = @(r,theta,phi) 0;
sigma_rtheta = @(r,theta,phi) 0;

D = @(theta,phi) [	sin(theta)^2*cos(phi)^2,        sin(theta)^2*sin(phi)^2,        cos(theta)^2,     	sin(2*theta)*sin(phi),      sin(2*theta)*cos(phi),  sin(theta)^2*sin(2*phi);
                    cos(theta)^2*cos(phi)^2,        cos(theta)^2*sin(phi)^2,        sin(theta)^2,     	-sin(2*theta)*sin(phi),     -sin(2*theta)*cos(phi), cos(theta)^2*sin(2*phi);
                    sin(phi)^2,                   	cos(phi)^2,                     0,               	0,                          0,                      -sin(2*phi);
                    -1/2*cos(theta)*sin(2*phi),     1/2*cos(theta)*sin(2*phi),      0,              	-sin(theta)*cos(phi),       sin(theta)*sin(phi),    cos(theta)*cos(2*phi);
                    -1/2*sin(theta)*sin(2*phi),     1/2*sin(theta)*sin(2*phi),      0,                	cos(theta)*cos(phi),        -cos(theta)*sin(phi),   sin(theta)*cos(2*phi);
                    1/2*sin(2*theta)*cos(phi)^2,	1/2*sin(2*theta)*sin(phi)^2,    -1/2*sin(2*theta),	cos(2*theta)*sin(phi),      cos(2*theta)*cos(phi),  1/2*sin(2*theta)*sin(2*phi)];
 
stress_u = @(r,theta,phi) D(theta,phi)\[sigma_rr(r,theta,phi);
                                        sigma_thetatheta(r,theta,phi);
                                        sigma_phiphi(r,theta,phi);
                                        sigma_thetaphi(r,theta,phi);
                                        sigma_rphi(r,theta,phi);
                                        sigma_rtheta(r,theta,phi)];
                                    
strain_u = @(r,theta,phi) C\stress_u(r,theta,phi);


stressMatrix = @(sigma) [sigma(1) sigma(6) sigma(5);
                         sigma(6) sigma(2) sigma(4);
                         sigma(5) sigma(4) sigma(3)];