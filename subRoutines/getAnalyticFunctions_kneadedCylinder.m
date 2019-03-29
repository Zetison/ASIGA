


u_analytic = @(r,theta,z) (1-2*nu)*(1+nu)/E*L/pi*[sin(c_1*pi*z/L); ...
                                                  sin(c_2*pi*z/L); ...
                                                  sin(c_3*pi*z/L)];

f = @(r,theta,z) [(1/2-nu)*c_1^2*pi*sin(c_1*pi*z/L)/L; ...
                  (1/2-nu)*c_2^2*pi*sin(c_2*pi*z/L)/L; ...
                  (1-nu)*c_3^2*sin(c_3*pi*z/L)*pi/L];

stress_u = @(x,y,z) [nu*c_3*cos(c_3*pi*z/L);
                         nu*c_3*cos(c_3*pi*z/L);
                         (1-nu)*c_3*cos(c_3*pi*z/L);
                         (1/2-nu)*c_2*cos(c_2*pi*z/L);
                         (1/2-nu)*c_1*cos(c_1*pi*z/L);
                         0];

stressMatrix = @(sigma) [sigma(1) sigma(6) sigma(5);
                         sigma(6) sigma(2) sigma(4);
                         sigma(5) sigma(4) sigma(3)];
                     
h_i = @(r,theta,z) -[nu*cos(c_3*pi*z/L)*c_3*cos(theta);
                     nu*cos(c_3*pi*z/L)*c_3*sin(theta);
                     0.5*(1-2*nu)*(cos(c_1*pi*z/L)*c_1*cos(theta)+cos(c_2*pi*z/L)*c_2*sin(theta))];
              
h_o = @(r,theta,z)  [nu*cos(c_3*pi*z/L)*c_3*cos(theta);
                     nu*cos(c_3*pi*z/L)*c_3*sin(theta);
                     0.5*(1-2*nu)*(cos(c_1*pi*z/L)*c_1*cos(theta)+cos(c_2*pi*z/L)*c_2*sin(theta))];
              

strain_u = @(radius,theta,z) [0;
                              0;
                              (1-2*nu)*(1+nu)*cos(c_3*pi*z/L)*c_3/E;
                              (1-2*nu)*(1+nu)*cos(c_2*pi*z/L)*c_2/E;
                              (1-2*nu)*(1+nu)*cos(c_1*pi*z/L)*c_1/E;
                              0];
                          
% f = @(r,theta,z) [-1+2*nu;
%                   -1+2*nu;
%                   -2+2*nu];
%                  
% h_i = @(r,theta,z) -[ -nu*(-2*z+L)*cos(theta);
%                   	  -nu*(-2*z+L)*sin(theta);
%                      0.5*(cos(theta)+sin(theta))*(-2*z+L)*(-1+2*nu)];
%               
% h_o = @(r,theta,z)  [ -nu*(-2*z+L)*cos(theta);
%                   	  -nu*(-2*z+L)*sin(theta);
%                      0.5*(cos(theta)+sin(theta))*(-2*z+L)*(-1+2*nu)];
%               
% 
% u_analytic = @(r,theta,z) (1-2*nu)*(1+nu)/E*[z*(z-L); ...
%                                              z*(z-L); ...
%                                              z*(z-L)];