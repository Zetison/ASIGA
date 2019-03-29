function [u, v, dudx, dudy, dudz, J] = numericalSolEval_newInfElementMethod(xi, eta, p_xi, p_eta, Xi, Eta, weights, ...
                                                    controlPts, U, Upsilon, k, N, r_o, noLayers, noDofs, element, D)
maxJacobian = 1e4;
if nargout <= 2
    R = NURBS2DBasis(xi, eta, p_xi, p_eta, Xi, Eta, weights);
else
    [R, dRdxi, dRdeta] = NURBS2DBasis(xi, eta,  p_xi, p_eta, Xi, Eta, weights);
end


uniqueXi = unique(Xi);
uniqueEta = unique(Eta);
noElementsXi = length(uniqueXi) - 1;
noElementsEta = length(uniqueEta) - 1;
xi_idx = findKnotSpan(noElementsXi, 0, xi, uniqueXi);
eta_idx = findKnotSpan(noElementsEta, 0, eta, uniqueEta);
e = xi_idx + noElementsXi*(eta_idx-1);
sctr = element(e,:);


vSurf = R*controlPts(sctr,:);

x = vSurf(1);
y = vSurf(2);
z = vSurf(3);

[chi, theta, phi] = evaluateProlateCoords(x,y,z,Upsilon);
r_arr = linspace(chi,r_o,noLayers);


u = zeros(1,noLayers);
for m = 1:N
    Q_m = zeros(1,noLayers);
    for mt = 1:N
        Q_m = Q_m + (chi./r_arr).^mt*D(m,mt);
    end
    phi_m = exp(1i*k*(r_arr-chi)).*Q_m;
    
    u = u + phi_m*(R*U(sctr + noDofs*(m-1)));    
end

v = [sqrt(r_arr.^2-Upsilon^2)*sin(theta)*cos(phi);
     sqrt(r_arr.^2-Upsilon^2)*sin(theta)*sin(phi);
     r_arr*cos(theta)];
 



