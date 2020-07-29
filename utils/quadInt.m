function I = quadInt(f,a,b)


[Q,W] = gaussQuad(64); 
I = 0;
detJ = 1/2*(b-a);
for i = 1:length(Q)
    pt = Q(i);
    x = a + (pt+1)*(b-a)/2;
    I = I + f(x)*W(i)*detJ;
end