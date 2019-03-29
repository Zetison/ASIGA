function [g,dg,S,dSdxi,dSdeta,D1,D2,D3] = getDepthRudderFunctions(b,l_ld,l_ud,b_ld,b_ud,h_d,delta_d,x_d,beta,s,c,part)

f_ld = @(xi) getNACA(xi,b_ld./l_ld,0);
f_ud = @(xi) getNACA(xi,b_ud./l_ud,0);
df_ld = @(xi) getNACA(xi,b_ld./l_ld,1);
df_ud = @(xi) getNACA(xi,b_ud/l_ud,1);


C_4 = c + b*cos(beta/2);
C_3 = b*sin(beta/2)-s;
C_2 = (2*C_4+C_3*tan(beta/2))/C_3^3;
C_1 = -(3*C_4+C_3*tan(beta/2))/C_3^2;

nPanels = 5;
dy = C_3/nPanels;
P = @(y) c + C_1*(y-s).^2 + C_2*(y-s).^3;


y_f = s+dy;

switch part
    case 'upper'
        pm = 1;
        Q1 = [0; s; c];
        Q2 = [0; y_f; c-b_ld/2];
        Q3 = [-1; y_f; c-b_ld/2];
    otherwise
        pm = -1;
        Q1 = [0; y_f+dy; P(y_f+dy)];
        Q2 = [0; y_f; c-b_ld/2];
        Q3 = [-1; y_f; c-b_ld/2];
end
normal = pm*cross(Q3-Q2,Q1-Q2);
D1 = normal(2);
D2 = normal(3);
D3 = D1*Q1(2) + D2*Q1(3);
g = @(xi) (D3 - D1*s - D2*(c - b_ld/2 + pm*l_ld*f_ld(xi)))./(D1*h_d + pm*D2*(l_ud.*f_ud(xi)-l_ld.*f_ld(xi)));
dg = @(xi) (-pm*D2*l_ld*df_ld(xi).*(D1*h_d + pm*D2*(l_ud.*f_ud(xi)-l_ld.*f_ld(xi))) - (D3 - D1*s - D2*(c - b_ld/2 + pm*l_ld*f_ld(xi))).*(pm*D2*(l_ud.*df_ud(xi)-l_ld.*df_ld(xi))))./(D1*h_d + pm*D2*(l_ud.*f_ud(xi)-l_ld.*f_ld(xi))).^2;
S = @(xi,eta) repmat([x_d;s;c-b_ld/2],1,size(xi,2)) + [-(l_ld.*xi.^2+eta.*(delta_d+(l_ud-l_ld).*xi.^2));
                                                      eta.*h_d;
                                                      pm*l_ld.*f_ld(xi) + pm*eta.*(l_ud.*f_ud(xi)-l_ld.*f_ld(xi))];

dSdxi = @(xi,eta) [-(2*l_ld.*xi+2*eta.*((l_ud-l_ld).*xi));
                  zeros(size(xi));
                  pm*l_ld.*df_ld(xi) + pm*eta.*(l_ud.*df_ud(xi)-l_ld.*df_ld(xi))];
dSdeta = @(xi,eta) [-((delta_d+(l_ud-l_ld).*xi.^2));
                      h_d*ones(size(xi));
                      pm*(l_ud.*f_ud(xi)-l_ld.*f_ld(xi))];

