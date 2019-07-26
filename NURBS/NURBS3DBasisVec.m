 function [R, dRdxi, dRdeta, dRdzeta] = NURBS3DBasisVec(xi, eta, zeta, ...
                                           p_xi, p_eta, p_zeta, Xi, Eta, Zeta, weights)
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

n_xi = length(Xi) - (p_xi+1);
n_eta = length(Eta) - (p_eta+1);
n_zeta = length(Zeta) - (p_zeta+1);

i1 = findKnotSpan(n_xi, p_xi, xi(1), Xi);
i2 = findKnotSpan(n_eta, p_eta, eta(1), Eta);
i3 = findKnotSpan(n_zeta, p_zeta, zeta(1), Zeta);

[N, dNdxi]   = Bspline_basisDers2(i1, xi, p_xi, Xi);
[M, dMdeta]  = Bspline_basisDers2(i2, eta, p_eta, Eta);
[L, dLdzeta] = Bspline_basisDers2(i3, zeta, p_zeta, Zeta);
    
noxi = numel(xi);

R = zeros(noxi, (p_xi+1)*(p_eta+1)*(p_zeta+1),class(xi));
if compDeriv
    dRdxi   = zeros(noxi, (p_xi+1)*(p_eta+1)*(p_zeta+1),class(xi));
    dRdeta  = zeros(noxi, (p_xi+1)*(p_eta+1)*(p_zeta+1),class(xi));
    dRdzeta = zeros(noxi, (p_xi+1)*(p_eta+1)*(p_zeta+1),class(xi));
end

W = zeros(noxi,1,class(xi));
if compDeriv
    dWdxi = zeros(noxi,1,class(xi));
    dWdeta = zeros(noxi,1,class(xi));
    dWdzeta = zeros(noxi,1,class(xi));
end

counter = 1;
for k3 = 1:p_zeta+1
    for k2 = 1:p_eta+1
        for k1 = 1:p_xi+1   
            weight = weights(counter);
            
            W = W + N(:,k1).*M(:,k2).*L(:,k3)*weight;
            if compDeriv
                dWdxi   = dWdxi   + dNdxi(:,k1).*M(:,k2)     .*L(:,k3)      .*weight;
                dWdeta  = dWdeta  + N(:,k1)    .*dMdeta(:,k2).*L(:,k3)      .*weight;
                dWdzeta = dWdzeta + N(:,k1)    .*M(:,k2)     .*dLdzeta(:,k3).*weight;
            end
            counter = counter + 1;
        end
    end
end

counter = 1;
for k3 = 1:p_zeta+1
    for k2 = 1:p_eta+1
        for k1 = 1:p_xi+1       
            fact = weights(counter)./(W.*W);
            
            NML = N(:,k1).*M(:,k2).*L(:,k3);
            R(:,counter) = NML.*fact.*W;
            if compDeriv
                dRdxi(:,counter)   = (dNdxi(:,k1)  .*M(:,k2).*L(:,k3).*W - NML.*dWdxi).*fact;
                dRdeta(:,counter)  = (dMdeta(:,k2) .*N(:,k1).*L(:,k3).*W - NML.*dWdeta).*fact;
                dRdzeta(:,counter) = (dLdzeta(:,k3).*N(:,k1).*M(:,k2).*W - NML.*dWdzeta).*fact;
            end
            counter = counter + 1;
        end
    end
end

