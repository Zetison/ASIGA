function p = exactKDT(k,P_inc,R)

p = 1i*P_inc./(4*k).*(1-(1+2*1i*k*R).*exp(-2*1i*k*R));

