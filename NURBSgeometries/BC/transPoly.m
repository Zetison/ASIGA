function P = transPoly(y,b,c,beta,s)
C_4 = c + b*cos(beta/2);
C_3 = b*sin(beta/2)-s;
C_2 = (2*C_4+C_3*tan(beta/2))/C_3^3;
C_1 = -(3*C_4+C_3*tan(beta/2))/C_3^2;

P = c + C_1*(y-s).^2 + C_2*(y-s).^3;