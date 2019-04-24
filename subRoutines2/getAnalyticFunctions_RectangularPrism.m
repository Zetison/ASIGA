
switch analyticFunction
    case 1
        u_analytic = @(x,y,z) [0;
                               x;
                               0];
                           
        stress_u = @(x,y,z) [0;
                             0;
                             0;
                             0;
                             0;
                             mu];
                                                   
        strain_u = @(x,y,z) [0;
                             0;
                             0;
                             0;
                             0;
                             1];
    case {2, 3}
        u_analytic = @(x,y,z) [(1/5)*x*y*(2*nu^2*x^2-6*nu^2*z^2-4*nu*x^2+15*nu*z^2-6*z^2);
                               (1/10)*x^2*(nu*x^2-3*nu*z^2-2*x^2+6*z^2);
                              -(3/5)*x^2*y*z*(nu-2)];
 
                             
        strain_u = @(x,y,z) [(1/5)*y*(2*nu^2*x^2-6*nu^2*z^2-4*nu*x^2+15*nu*z^2-6*z^2)+(1/5)*x*y*(4*nu^2*x-8*nu*x);
                             0;
                            -(3/5)*x^2*y*(nu-2);
                            (1/10)*x^2*(-6*nu*z+12*z)-(3/5)*x^2*z*(nu-2);
                            -(6/5)*x*y*z*(nu-2)+(1/5)*x*y*(-12*nu^2*z+30*nu*z-12*z);
                             (1/5)*x*(2*nu^2*x^2-6*nu^2*z^2-4*nu*x^2+15*nu*z^2-6*z^2)+(1/5)*x*(nu*x^2-3*nu*z^2-2*x^2+6*z^2)+(1/10)*x^2*(2*nu*x-4*x)];
                           
        stress_u = @(x,y,z) [(3/5)*(nu^2*x^2-nu^2*z^2-2*nu*x^2+3*nu*z^2-2*z^2)*E*y/(1+nu);
                                 -(3/5)*(nu*x^2-nu*z^2-2*x^2+2*z^2)*E*nu*y/(1+nu);
                                 -(3/5)*(nu^2*x^2-nu^2*z^2-nu*x^2+2*nu*z^2-2*x^2)*E*y/(1+nu);
                                 -(3/5)*(nu-2)*z*x^2*E/(1+nu);
                                 -(6/5)*(nu-2)*nu*z*y*x*E/(1+nu);
                                 (1/5)*(nu^2*x^2-3*nu^2*z^2-nu*x^2+6*nu*z^2-2*x^2)*x*E/(1+nu)];
                      
    case 4
        u_analytic = @(x,y,z) [x*(x-wx)*y*(y-wy)*z*(z-wz);
                               x*(x-wx)*y*(y-wy)*z*(z-wz);
                               x*(x-wx)*y*(y-wy)*z*(z-wz)];

        strain_u = @(x,y,z) [(x-wx)*y*(y-wy)*z*(z-wz)+x*y*(y-wy)*z*(z-wz);
                    x*(x-wx)*(y-wy)*z*(z-wz)+x*(x-wx)*y*z*(z-wz);
                    x*(x-wx)*y*(y-wy)*(z-wz)+x*(x-wx)*y*(y-wy)*z;
                    x*(x-wx)*y*(y-wy)*(z-wz)+x*(x-wx)*y*(y-wy)*z+x*(x-wx)*(y-wy)*z*(z-wz)+x*(x-wx)*y*z*(z-wz);
                    (x-wx)*y*(y-wy)*z*(z-wz)+x*y*(y-wy)*z*(z-wz)+x*(x-wx)*y*(y-wy)*(z-wz)+x*(x-wx)*y*(y-wy)*z;
                    x*(x-wx)*(y-wy)*z*(z-wz)+x*(x-wx)*y*z*(z-wz)+(x-wx)*y*(y-wy)*z*(z-wz)+x*y*(y-wy)*z*(z-wz)];

        f = @(x,y,z) [-1/2*1/(1+nu)/(-1+2*nu)*(E*(4*nu*wx*wy*x*y+4*nu*wx*wz*x*z-4*nu*wx*x*y^2-4*nu*wx*x*z^2+4*nu*wy*wz*y*z-4*nu*wy*x^2*y-4*nu*wy*y*z^2-4*nu*wz*x^2*z-4*nu*wz*y^2*z+4*nu*x^2*y^2+4*nu*x^2*z^2+4*nu*y^2*z^2+wx*wy*wz*y+wx*wy*wz*z-2*wx*wy*x*y-2*wx*wy*y*z-wx*wy*z^2-2*wx*wz*x*z-wx*wz*y^2-2*wx*wz*y*z+2*wx*x*y^2+2*wx*x*z^2+2*wx*y^2*z+2*wx*y*z^2-2*wy*wz*x*y-2*wy*wz*x*z-4*wy*wz*y*z+2*wy*x^2*y+4*wy*x*y*z+2*wy*x*z^2+4*wy*y*z^2+2*wz*x^2*z+2*wz*x*y^2+4*wz*x*y*z+4*wz*y^2*z-2*x^2*y^2-2*x^2*z^2-4*x*y^2*z-4*x*y*z^2-4*y^2*z^2));
                      -1/2*1/(1+nu)/(-1+2*nu)*(E*(4*nu*wx*wy*x*y+4*nu*wx*wz*x*z-4*nu*wx*x*y^2-4*nu*wx*x*z^2+4*nu*wy*wz*y*z-4*nu*wy*x^2*y-4*nu*wy*y*z^2-4*nu*wz*x^2*z-4*nu*wz*y^2*z+4*nu*x^2*y^2+4*nu*x^2*z^2+4*nu*y^2*z^2+wx*wy*wz*x+wx*wy*wz*z-2*wx*wy*x*y-2*wx*wy*x*z-wx*wy*z^2-2*wx*wz*x*y-4*wx*wz*x*z-2*wx*wz*y*z+2*wx*x*y^2+4*wx*x*y*z+4*wx*x*z^2+2*wx*y*z^2-wy*wz*x^2-2*wy*wz*x*z-2*wy*wz*y*z+2*wy*x^2*y+2*wy*x^2*z+2*wy*x*z^2+2*wy*y*z^2+2*wz*x^2*y+4*wz*x^2*z+4*wz*x*y*z+2*wz*y^2*z-2*x^2*y^2-4*x^2*y*z-4*x^2*z^2-4*x*y*z^2-2*y^2*z^2));
                      -1/2*1/(1+nu)/(-1+2*nu)*(E*(4*nu*wx*wy*x*y+4*nu*wx*wz*x*z-4*nu*wx*x*y^2-4*nu*wx*x*z^2+4*nu*wy*wz*y*z-4*nu*wy*x^2*y-4*nu*wy*y*z^2-4*nu*wz*x^2*z-4*nu*wz*y^2*z+4*nu*x^2*y^2+4*nu*x^2*z^2+4*nu*y^2*z^2+wx*wy*wz*x+wx*wy*wz*y-4*wx*wy*x*y-2*wx*wy*x*z-2*wx*wy*y*z-2*wx*wz*x*y-2*wx*wz*x*z-wx*wz*y^2+4*wx*x*y^2+4*wx*x*y*z+2*wx*x*z^2+2*wx*y^2*z-wy*wz*x^2-2*wy*wz*x*y-2*wy*wz*y*z+4*wy*x^2*y+2*wy*x^2*z+4*wy*x*y*z+2*wy*y*z^2+2*wz*x^2*y+2*wz*x^2*z+2*wz*x*y^2+2*wz*y^2*z-4*x^2*y^2-4*x^2*y*z-2*x^2*z^2-4*x*y^2*z-2*y^2*z^2))];

                
    case 5
        u_analytic = @(x,y,z) [x^2*(x-wx)*y*(y-wy)*z*(z-wz);
                               x^2*(x-wx)*y*(y-wy)*z*(z-wz);
                               x^2*(x-wx)*y*(y-wy)*z*(z-wz)];

        strain_u = @(x,y,z) [2*x*(x-wx)*y*(y-wy)*z*(z-wz)+x^2*y*(y-wy)*z*(z-wz);
                             x^2*(x-wx)*(y-wy)*z*(z-wz)+x^2*(x-wx)*y*z*(z-wz);
                             x^2*(x-wx)*y*(y-wy)*(z-wz)+x^2*(x-wx)*y*(y-wy)*z;
                             x^2*(x-wx)*y*(y-wy)*(z-wz)+x^2*(x-wx)*y*(y-wy)*z+x^2*(x-wx)*(y-wy)*z*(z-wz)+x^2*(x-wx)*y*z*(z-wz);
                             2*x*(x-wx)*y*(y-wy)*z*(z-wz)+x^2*y*(y-wy)*z*(z-wz)+x^2*(x-wx)*y*(y-wy)*(z-wz)+x^2*(x-wx)*y*(y-wy)*z;
                             x^2*(x-wx)*(y-wy)*z*(z-wz)+x^2*(x-wx)*y*z*(z-wz)+2*x*(x-wx)*y*(y-wy)*z*(z-wz)+x^2*y*(y-wy)*z*(z-wz)];

        f = @(x,y,z) [1/2*1/(1+nu)/(-1+2*nu)*(E*(4*nu*wx*wy*wz*y*z-4*nu*wx*wy*x^2*y-4*nu*wx*wy*y*z^2-4*nu*wx*wz*x^2*z-4*nu*wx*wz*y^2*z+4*nu*wx*x^2*y^2+4*nu*wx*x^2*z^2+4*nu*wx*y^2*z^2-12*nu*wy*wz*x*y*z+4*nu*wy*x^3*y+12*nu*wy*x*y*z^2+4*nu*wz*x^3*z+12*nu*wz*x*y^2*z-4*nu*x^3*y^2-4*nu*x^3*z^2-12*nu*x*y^2*z^2-2*wx*wy*wz*x*y-2*wx*wy*wz*x*z-4*wx*wy*wz*y*z+2*wx*wy*x^2*y+4*wx*wy*x*y*z+2*wx*wy*x*z^2+4*wx*wy*y*z^2+2*wx*wz*x^2*z+2*wx*wz*x*y^2+4*wx*wz*x*y*z+4*wx*wz*y^2*z-2*wx*x^2*y^2-2*wx*x^2*z^2-4*wx*x*y^2*z-4*wx*x*y*z^2-4*wx*y^2*z^2+3*wy*wz*x^2*y+3*wy*wz*x^2*z+12*wy*wz*x*y*z-2*wy*x^3*y-6*wy*x^2*y*z-3*wy*x^2*z^2-12*wy*x*y*z^2-2*wz*x^3*z-3*wz*x^2*y^2-6*wz*x^2*y*z-12*wz*x*y^2*z+2*x^3*y^2+2*x^3*z^2+6*x^2*y^2*z+6*x^2*y*z^2+12*x*y^2*z^2));
                      1/2*1/(1+nu)/(-1+2*nu)*(E*(4*nu*wx*wy*wz*y*z-4*nu*wx*wy*x^2*y-4*nu*wx*wy*y*z^2-4*nu*wx*wz*x^2*z-4*nu*wx*wz*y^2*z+4*nu*wx*x^2*y^2+4*nu*wx*x^2*z^2+4*nu*wx*y^2*z^2-12*nu*wy*wz*x*y*z+4*nu*wy*x^3*y+12*nu*wy*x*y*z^2+4*nu*wz*x^3*z+12*nu*wz*x*y^2*z-4*nu*x^3*y^2-4*nu*x^3*z^2-12*nu*x*y^2*z^2-wx*wy*wz*x^2-2*wx*wy*wz*x*z-2*wx*wy*wz*y*z+2*wx*wy*x^2*y+2*wx*wy*x^2*z+2*wx*wy*x*z^2+2*wx*wy*y*z^2+2*wx*wz*x^2*y+4*wx*wz*x^2*z+4*wx*wz*x*y*z+2*wx*wz*y^2*z-2*wx*x^2*y^2-4*wx*x^2*y*z-4*wx*x^2*z^2-4*wx*x*y*z^2-2*wx*y^2*z^2+wy*wz*x^3+3*wy*wz*x^2*z+6*wy*wz*x*y*z-2*wy*x^3*y-2*wy*x^3*z-3*wy*x^2*z^2-6*wy*x*y*z^2-2*wz*x^3*y-4*wz*x^3*z-6*wz*x^2*y*z-6*wz*x*y^2*z+2*x^3*y^2+4*x^3*y*z+4*x^3*z^2+6*x^2*y*z^2+6*x*y^2*z^2));
                      1/2*1/(1+nu)/(-1+2*nu)*(E*(4*nu*wx*wy*wz*y*z-4*nu*wx*wy*x^2*y-4*nu*wx*wy*y*z^2-4*nu*wx*wz*x^2*z-4*nu*wx*wz*y^2*z+4*nu*wx*x^2*y^2+4*nu*wx*x^2*z^2+4*nu*wx*y^2*z^2-12*nu*wy*wz*x*y*z+4*nu*wy*x^3*y+12*nu*wy*x*y*z^2+4*nu*wz*x^3*z+12*nu*wz*x*y^2*z-4*nu*x^3*y^2-4*nu*x^3*z^2-12*nu*x*y^2*z^2-wx*wy*wz*x^2-2*wx*wy*wz*x*y-2*wx*wy*wz*y*z+4*wx*wy*x^2*y+2*wx*wy*x^2*z+4*wx*wy*x*y*z+2*wx*wy*y*z^2+2*wx*wz*x^2*y+2*wx*wz*x^2*z+2*wx*wz*x*y^2+2*wx*wz*y^2*z-4*wx*x^2*y^2-4*wx*x^2*y*z-2*wx*x^2*z^2-4*wx*x*y^2*z-2*wx*y^2*z^2+wy*wz*x^3+3*wy*wz*x^2*y+6*wy*wz*x*y*z-4*wy*x^3*y-2*wy*x^3*z-6*wy*x^2*y*z-6*wy*x*y*z^2-2*wz*x^3*y-2*wz*x^3*z-3*wz*x^2*y^2-6*wz*x*y^2*z+4*x^3*y^2+4*x^3*y*z+2*x^3*z^2+6*x^2*y^2*z+6*x*y^2*z^2))];

           
    case 6
        u_analytic = @(x,y,z) [sin(4*pi*x/wx)*sin(3*pi*y/wy)*sin(2*pi*z/wz);
                               sin(4*pi*x/wx)*sin(3*pi*y/wy)*sin(2*pi*z/wz);
                               sin(4*pi*x/wx)*sin(3*pi*y/wy)*sin(2*pi*z/wz)];
                
        strain_u = @(x,y,z) [4*pi*cos(4*pi*x/wx)*sin(3*pi*y/wy)*sin(2*pi*z/wz)/wx;
                             3*sin(4*pi*x/wx)*pi*cos(3*pi*y/wy)*sin(2*pi*z/wz)/wy;
                             2*sin(4*pi*x/wx)*sin(3*pi*y/wy)*pi*cos(2*pi*z/wz)/wz;
                             2*sin(4*pi*x/wx)*sin(3*pi*y/wy)*pi*cos(2*pi*z/wz)/wz+3*sin(4*pi*x/wx)*pi*cos(3*pi*y/wy)*sin(2*pi*z/wz)/wy;
                             4*pi*cos(4*pi*x/wx)*sin(3*pi*y/wy)*sin(2*pi*z/wz)/wx+2*sin(4*pi*x/wx)*sin(3*pi*y/wy)*pi*cos(2*pi*z/wz)/wz;
                             3*sin(4*pi*x/wx)*pi*cos(3*pi*y/wy)*sin(2*pi*z/wz)/wy+4*pi*cos(4*pi*x/wx)*sin(3*pi*y/wy)*sin(2*pi*z/wz)/wx];

        f = @(x,y,z) [1/2*1/(1+nu)/(-1+2*nu)/wx^2/wy^2/wz^2*(E*pi^2*(8*sin(4*pi*x/wx)*sin(3*pi*y/wy)*sin(2*pi*z/wz)*wx^2*wy^2*nu+18*sin(4*pi*x/wx)*sin(3*pi*y/wy)*sin(2*pi*z/wz)*nu*wx^2*wz^2+32*sin(4*pi*x/wx)*sin(3*pi*y/wy)*sin(2*pi*z/wz)*nu*wy^2*wz^2-4*sin(4*pi*x/wx)*sin(3*pi*y/wy)*sin(2*pi*z/wz)*wx^2*wy^2-9*sin(4*pi*x/wx)*sin(3*pi*y/wy)*sin(2*pi*z/wz)*wx^2*wz^2-32*sin(4*pi*x/wx)*sin(3*pi*y/wy)*sin(2*pi*z/wz)*wy^2*wz^2+8*sin(3*pi*y/wy)*cos(4*pi*x/wx)*cos(2*pi*z/wz)*wx*wy^2*wz+12*sin(2*pi*z/wz)*cos(4*pi*x/wx)*cos(3*pi*y/wy)*wx*wy*wz^2)); 
                      1/(1+nu)/(-1+2*nu)/wx^2/wy^2/wz^2*(E*pi^2*(4*sin(4*pi*x/wx)*sin(3*pi*y/wy)*sin(2*pi*z/wz)*wx^2*wy^2*nu+9*sin(4*pi*x/wx)*sin(3*pi*y/wy)*sin(2*pi*z/wz)*nu*wx^2*wz^2+16*sin(4*pi*x/wx)*sin(3*pi*y/wy)*sin(2*pi*z/wz)*nu*wy^2*wz^2-2*sin(4*pi*x/wx)*sin(3*pi*y/wy)*sin(2*pi*z/wz)*wx^2*wy^2-9*sin(4*pi*x/wx)*sin(3*pi*y/wy)*sin(2*pi*z/wz)*wx^2*wz^2-8*sin(4*pi*x/wx)*sin(3*pi*y/wy)*sin(2*pi*z/wz)*wy^2*wz^2+3*sin(4*pi*x/wx)*cos(2*pi*z/wz)*cos(3*pi*y/wy)*wx^2*wy*wz+6*sin(2*pi*z/wz)*cos(4*pi*x/wx)*cos(3*pi*y/wy)*wx*wy*wz^2));
                      1/2*1/(1+nu)/(-1+2*nu)/wx^2/wy^2/wz^2*(E*pi^2*(8*sin(4*pi*x/wx)*sin(3*pi*y/wy)*sin(2*pi*z/wz)*wx^2*wy^2*nu+18*sin(4*pi*x/wx)*sin(3*pi*y/wy)*sin(2*pi*z/wz)*nu*wx^2*wz^2+32*sin(4*pi*x/wx)*sin(3*pi*y/wy)*sin(2*pi*z/wz)*nu*wy^2*wz^2-8*sin(4*pi*x/wx)*sin(3*pi*y/wy)*sin(2*pi*z/wz)*wx^2*wy^2-9*sin(4*pi*x/wx)*sin(3*pi*y/wy)*sin(2*pi*z/wz)*wx^2*wz^2-16*sin(4*pi*x/wx)*sin(3*pi*y/wy)*sin(2*pi*z/wz)*wy^2*wz^2+6*sin(4*pi*x/wx)*cos(2*pi*z/wz)*cos(3*pi*y/wy)*wx^2*wy*wz+8*sin(3*pi*y/wy)*cos(4*pi*x/wx)*cos(2*pi*z/wz)*wx*wy^2*wz))];
    


                
end

stressMatrix = @(sigma) [sigma(1) sigma(6) sigma(5);
                         sigma(6) sigma(2) sigma(4);
                         sigma(5) sigma(4) sigma(3)];