function [g,dg,S,dSdxi,dSdeta] = getMainRudderFunctions(b,g2,L,l_lm,l_um,b_lm,b_um,h_m,delta_m,alpha,x_m)

h = g2*tan(alpha./2); % = g2./tan((pi-alpha)./2)

x_c = -(L+g2+(b-h).*cot(alpha));

f_lm = @(xi) getNACA(xi,b_lm/l_lm,0);
f_um = @(xi) getNACA(xi,b_um/l_um,0);
df_lm = @(xi) getNACA(xi,b_lm/l_lm,1);
df_um = @(xi) getNACA(xi,b_um/l_um,1);

C_a = @(xi) (l_um*f_um(xi) - l_lm*f_lm(xi)).^2+h_m^2 - delta_m^2*(1-xi.^2).^2*tan(alpha)^2;
C_b = @(xi) 2*l_lm*f_lm(xi).*(l_um*f_um(xi)-l_lm*f_lm(xi)) + 2*tan(alpha)^2*(x_m-l_lm*xi.^2-x_c)*delta_m.*(1-xi.^2);
C_c = @(xi) (l_lm*f_lm(xi)).^2 - tan(alpha)^2*(x_m-l_lm*xi.^2-x_c).^2;
g = @(xi) (-C_b(xi)+sqrt(C_b(xi).^2-4*C_a(xi).*C_c(xi)))./(2*C_a(xi));

dC_a = @(xi) 2*(l_um*f_um(xi) - l_lm*f_lm(xi)).*(l_um*df_um(xi) - l_lm*df_lm(xi)) + delta_m^2*4*xi.*(1-xi.^2)*tan(alpha)^2;
dC_b = @(xi) 2*l_lm*df_lm(xi).*(l_um*f_um(xi)-l_lm*f_lm(xi)) + 2*l_lm*f_lm(xi).*(l_um*df_um(xi)-l_lm*df_lm(xi)) - 4*tan(alpha)^2*l_lm*xi*delta_m.*(1-xi.^2) - 4*tan(alpha)^2*(x_m-l_lm*xi.^2-x_c)*delta_m.*xi;
dC_c = @(xi) 2*l_lm^2*f_lm(xi).*df_lm(xi) + 4*tan(alpha)^2*(x_m-l_lm*xi.^2-x_c)*l_lm.*xi;


dg = @(xi) ((-dC_b(xi) + 1./sqrt(C_b(xi).^2-4*C_a(xi).*C_c(xi)).*(C_b(xi).*dC_b(xi)-2*(dC_a(xi).*C_c(xi) + C_a(xi).*dC_c(xi))))*2.*C_a(xi) ...
            - (-C_b(xi) + sqrt(C_b(xi).^2-4*C_a(xi).*C_c(xi)))*2.*dC_a(xi))./(2*C_a(xi)).^2;

S = @(xi,eta) repmat([x_m,0,0],size(xi,1),1) + [-(l_lm.*xi.^2+eta*delta_m.*(1-xi.^2)), ...
                                                  l_lm.*f_lm(xi) + eta.*(l_um.*f_um(xi)-l_lm.*f_lm(xi)), ...
                                                  eta.*h_m];

dSdxi = @(xi,eta) [-(2*l_lm.*xi-2*eta*delta_m.*xi), ...
                   l_lm.*df_lm(xi) + eta.*(l_um.*df_um(xi)-l_lm.*df_lm(xi)), ...
                   zeros(size(xi,1),1)];
dSdeta = @(xi,eta) [-delta_m.*(1-xi.^2), ...
                      l_um.*f_um(xi)-l_lm.*f_lm(xi), ...
                      h_m*ones(size(xi,1),1)];