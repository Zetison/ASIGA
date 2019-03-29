
% u_analytic = @(r,theta,z) 1/E*p_i*R_i^2/(R_o^2-R_i^2)*((1+nu)*(1-2*nu)*r + R_o^2*(1+nu)/r);
u_analytic = @(r,theta,z) (1/2)*((R_i^2*p_i-R_o^2*p_o)*r/(lambda+mu)+R_i^2*R_o^2*(p_i-p_o)/(mu*r))/(-R_i^2+R_o^2);

sigma_rr = @(r,theta,z) (R_i^2*p_i - R_o^2*p_o)/(R_o^2-R_i^2) + (R_i^2*R_o^2*p_o - R_i^2*R_o^2*p_i)/(R_o^2-R_i^2)/r^2;
sigma_thetatheta = @(r,theta,z) (R_i^2*p_i - R_o^2*p_o)/(R_o^2-R_i^2) + (R_i^2*R_o^2*p_i - R_i^2*R_o^2*p_o)/(R_o^2-R_i^2)/r^2;
sigma_zz = @(r,theta,z) 2*nu*(R_i^2*p_i-R_o^2*p_o)/(R_o^2-R_i^2);
sigma_thetaz = @(r,theta,z) 0;
sigma_rz = @(r,theta,z) 0;
sigma_rtheta = @(r,theta,z) 0;


stress_u = @(r,theta,z) [(cos(theta))^2*sigma_rr(r,theta,z) + (sin(theta))^2*sigma_thetatheta(r,theta,z) - sin(2*theta)*sigma_rtheta(r,theta,z);
                              (sin(theta))^2*sigma_rr(r,theta,z) + (cos(theta))^2*sigma_thetatheta(r,theta,z) + sin(2*theta)*sigma_rtheta(r,theta,z);
                              sigma_zz(r,theta,z);
                              sigma_thetaz(r,theta,z)*cos(theta) + sigma_rz(r,theta,z)*sin(theta);
                              -sigma_thetaz(r,theta,z)*sin(theta) + sigma_rz(r,theta,z)*cos(theta);
                              1/2*sin(2*theta)*sigma_rr(r,theta,z) - 1/2*sin(2*theta)*sigma_thetatheta(r,theta,z) - cos(2*theta)*sigma_rtheta(r,theta,z)];
strain_u = @(r,theta,z) C\stress_u(r,theta,z);


stressMatrix = @(sigma) [sigma(1) sigma(6) sigma(5);
                         sigma(6) sigma(2) sigma(4);
                         sigma(5) sigma(4) sigma(3)];