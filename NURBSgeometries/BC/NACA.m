function [f,dfxi,dfdxi2] = NACA(b,l)

f = @(xi) getNACA(xi,b/l,0);
dfxi = @(xi) getNACA(xi,b/l,1);
dfdxi2 = @(xi) getNACA(xi,b/l,2);