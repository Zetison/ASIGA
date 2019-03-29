function [R, dRdxi, dRdeta] = NURBS2DBasisVec(xi, eta, p_xi, p_eta, Xi, Eta, weights)

n_xi = length(Xi) - (p_xi+1);
n_eta = length(Eta) - (p_eta+1);

i1 = findKnotSpan(n_xi, p_xi, xi(1), Xi);
i2 = findKnotSpan(n_eta, p_eta, eta(2), Eta);

[N, dNdxi] = Bspline_basisDers2(i1, xi, p_xi, Xi);
[M, dMdeta] = Bspline_basisDers2(i2, eta, p_eta, Eta);

noxi = numel(xi);

R = zeros(noxi, (p_xi+1)*(p_eta+1));
dRdxi = zeros(noxi, (p_xi+1)*(p_eta+1));
dRdeta = zeros(noxi, (p_xi+1)*(p_eta+1));

W = zeros(noxi,1);
dWdxi = zeros(noxi,1);
dWdeta = zeros(noxi,1);

counter = 1;
for k2 = 1:p_eta+1
    for k1 = 1:p_xi+1    
        weight = weights(counter);

        W       = W       + N(:,k1)    .*M(:,k2)     *weight;
        dWdxi   = dWdxi   + dNdxi(:,k1).*M(:,k2)     *weight;
        dWdeta  = dWdeta  + N(:,k1)    .*dMdeta(:,k2)*weight;
        counter = counter + 1;
    end
end

counter = 1;
for k2 = 1:p_eta+1
    for k1 = 1:p_xi+1     
        fact = weights(counter)./(W.*W);

        NM = N(:,k1).*M(:,k2);
        R(:,counter) = NM.*fact.*W;

        dRdxi(:,counter)   = (dNdxi(:,k1)  .*M(:,k2).*W - NM.*dWdxi).*fact;
        dRdeta(:,counter)  = (dMdeta(:,k2) .*N(:,k1).*W - NM.*dWdeta).*fact;
        counter = counter + 1;
    end
end

