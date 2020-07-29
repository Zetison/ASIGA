 function [R, dRdxi, dRdeta, dRdzeta] = NURBS3DBasis(xi, eta, zeta, ...
                                           p, q, r, Xi, Eta, Zeta, weights)
error('Depricated. Use NURBSbasis instead')
% This routine compute the (p+1)(q+1)(r+1) nonzero NURBS functions 
% and corresponding derivatives at (xi, eta, zeta)

% Input
%       (xi,eta,zeta):      evaluation point
%       p,q,r:              NURBS degrees
%       Xi, Eta, Zeta:      knot vectors
%       weights:            NURBS weights

% Output
%       R, dRdxi, dRdeta, dRdzeta:      array of the (p+1)(q+1)(r+1) NURBS 
%                                       functions and its derivativeswhich 
%                                       are nonzero at (xi,eta,zeta)
if nargout == 1
    compDeriv = false;
else
    compDeriv = true;
end

n = length(Xi) - (p+1);
m = length(Eta) - (q+1);
l = length(Zeta) - (r+1);

i1 = findKnotSpan(n, p, xi, Xi);
i2 = findKnotSpan(m, q, eta, Eta);
i3 = findKnotSpan(l, r, zeta, Zeta);

[N, dNdxi]   = Bspline_basis(i1, xi, p, Xi, compDeriv);
[M, dMdeta]  = Bspline_basis(i2, eta, q, Eta, compDeriv);
[L, dLdzeta] = Bspline_basis(i3, zeta, r, Zeta, compDeriv);
    

R = zeros(1, (p+1)*(q+1)*(r+1),class(xi));
if compDeriv
    dRdxi   = zeros(1, (p+1)*(q+1)*(r+1),class(xi));
    dRdeta  = zeros(1, (p+1)*(q+1)*(r+1),class(xi));
    dRdzeta = zeros(1, (p+1)*(q+1)*(r+1),class(xi));
end

W = 0;
if compDeriv
    dWdxi   = 0;
    dWdeta  = 0;
    dWdzeta = 0;
end

counter = 1;
for k3 = 1:r+1
    for k2 = 1:q+1
        for k1 = 1:p+1   
            weight = weights(counter);
            
            W = W + N(k1)*M(k2)*L(k3)*weight;
            if compDeriv
                dWdxi   = dWdxi   + dNdxi(k1)*M(k2)     *L(k3)      *weight;
                dWdeta  = dWdeta  + N(k1)    *dMdeta(k2)*L(k3)      *weight;
                dWdzeta = dWdzeta + N(k1)    *M(k2)     *dLdzeta(k3)*weight;
            end
            counter = counter + 1;
        end
    end
end

counter = 1;
for k3 = 1:r+1
    for k2 = 1:q+1
        for k1 = 1:p+1       
            fact = weights(counter)/(W*W);
            
            NML = N(k1)*M(k2)*L(k3);
            R(counter) = NML*fact*W;
            if compDeriv
                dRdxi(counter)   = (dNdxi(k1)  *M(k2)*L(k3)*W - NML*dWdxi)*fact;
                dRdeta(counter)  = (dMdeta(k2) *N(k1)*L(k3)*W - NML*dWdeta)*fact;
                dRdzeta(counter) = (dLdzeta(k3)*N(k1)*M(k2)*W - NML*dWdzeta)*fact;
            end
            counter = counter + 1;
        end
    end
end

