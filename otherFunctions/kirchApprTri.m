function p = kirchApprTri(tri,P,X,varCol)

P_inc = varCol.P_inc;
k = varCol.k;
Eps = 1e2*eps;
P1 = P(tri(:,1),:);
P2 = P(tri(:,2),:);
P3 = P(tri(:,3),:);

D0 = -2*1i*(P3*X.')*k;
D1 = -2*1i*((P1-P3)*X.')*k;
D2 = -2*1i*((P2-P3)*X.')*k;

g = (D1.*(1-exp(D2))+D2.*(exp(D1)-1))./D1./D2./(D1-D2);

indices1 = abs(D1) < Eps;
temp = D2(indices1);
g(indices1) = -(1+temp-exp(temp))./temp.^2;
indices2 = abs(D2) < Eps;
temp = D1(indices2);
g(indices2) = -(1+temp-exp(temp))./temp.^2;
indices3 = abs(D1-D2) < Eps;
temp = D2(indices3);
g(indices3) = (1-exp(temp)+temp.*exp(temp))./temp.^2;
g(and(indices1,indices2)) = 0.5;
normals = cross(P2-P1,P3-P1,2);
areas = norm2(normals); % (multiplied by 2)
n = normals./areas(:,[1,1,1]);
nX = n*X.';
areas(nX < 0) = 0;
p = -1i*k*P_inc/(2*pi).*sum(areas.*nX.*g.*exp(D0),1);

