
% k_vec = [5 4 3];
% k_wn = norm(k_vec);
% Grad = @(x,y,z) conj([1i*k_vec(1);
%                  1i*k_vec(2);
%                  1i*k_vec(3)]*exp(1i*dot(k_vec,[x,y,z])));
% 
% u_analytic = @(x,y,z) exp(1i*(k_vec(1)*x+k_vec(2)*y+k_vec(3)*z));

k_wn = 1;
c_f = 1500;
omega = k_wn*c_f;
             
             
Acoeff = [1, 1i+1];
Bcoeff = [1-1i, 1i];
Ccoeff = [0.5, 1i+1];
Dcoeff = [1-1i, 1i];

Ecoeff = [1-1i, 1i;
          1,    1];
 
Fcoeff = [1+1i, 1;
          1i,    -1];
      
u_analytic = @(x,y,z) general3DSolutionHelmholtz(x,y,z,k_wn, Acoeff, Bcoeff, Ccoeff, Dcoeff, Ecoeff, Fcoeff);
Grad =  @(x,y,z) general3DSolutionHelmholtz_grad(x,y,z,k_wn, Acoeff, Bcoeff, Ccoeff, Dcoeff, Ecoeff, Fcoeff);


