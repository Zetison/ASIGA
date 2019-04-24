
% u_analytic = @(r,theta,phi) 1/E*p_i*R_i^2/(R_o^2-R_i^2)*((1+nu)*(1-2*nu)*r + R_o^2*(1+nu)/r);
u_analytic = @(r,theta,phi) 1/(R_o^3-R_i^3)*((R_i^3*p_i - R_o^3*p_o)*r/(3*lambda+2*mu) + 1/4*R_i^3*R_o^3*(p_i-p_o)/(mu*r^2));

sigma_rr = @(r,theta,phi) 1/(R_o^3-R_i^3)*(R_i^3*R_o^3*(p_o-p_i)/r^3 + R_i^3*p_i - R_o^3*p_o);
sigma_thetatheta = @(r,theta,phi) 1/2*1/(R_o^3-R_i^3)*(R_i^3*R_o^3*(p_i-p_o)/r^3 + 2*R_i^3*p_i - 2*R_o^3*p_o);
sigma_phiphi = @(r,theta,phi) 1/2*1/(R_o^3-R_i^3)*(R_i^3*R_o^3*(p_i-p_o)/r^3 + 2*R_i^3*p_i - 2*R_o^3*p_o);
sigma_thetaphi = @(r,theta,phi) 0;
sigma_rphi = @(r,theta,phi) 0;
sigma_rtheta = @(r,theta,phi) 0;

D = @(theta,phi) [	sin(theta)^2*cos(phi)^2,       sin(theta)^2*sin(phi)^2,        cos(theta)^2,     	sin(2*theta)*sin(phi),      sin(2*theta)*cos(phi),  sin(theta)^2*sin(2*phi);
                    cos(theta)^2*cos(phi)^2,       cos(theta)^2*sin(phi)^2,        sin(theta)^2,     	-sin(2*theta)*sin(phi),     -sin(2*theta)*cos(phi), cos(theta)^2*sin(2*phi);
                    sin(phi)^2,                    cos(phi)^2,                     0,               	0,                          0,                      -sin(2*phi);
                    -1/2*cos(theta)*sin(2*phi),	1/2*cos(theta)*sin(2*phi),      0,              	-sin(theta)*cos(phi),       sin(theta)*sin(phi),    cos(theta)*cos(2*phi);
                    -1/2*sin(theta)*sin(2*phi),	1/2*sin(theta)*sin(2*phi),      0,                	cos(theta)*cos(phi),        -cos(theta)*sin(phi),   sin(theta)*cos(2*phi);
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