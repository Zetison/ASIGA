function SS = getBeTSSiSail(ss,eta,a,c,l_ls,l_us,b_ls,b_us,h_s,delta_s,s2)
t = b_ls/l_ls;
xi = invertNACA(ss,s2,t);

T = [a-19,0,c].';
[f_u,dfxi_u,dfdxi2_u] = NACA(b_ls,l_ls);
[f_o,dfxi_o,dfdxi2_o] = NACA(b_us,l_us);
X = rudder(l_ls,l_us,h_s,delta_s,eye(3),T,1,f_u,dfxi_u,dfdxi2_u,f_o,dfxi_o,dfdxi2_o);
SS = X([xi;eta]);
        
        
        
        
        
        
        