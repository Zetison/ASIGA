function calculateSolutionAlongRay(varCol, U, k, xi, eta)
% close all
Xi = varCol.nurbs.knots{1};
Eta = varCol.nurbs.knots{2};
Zeta = varCol.nurbs.knots{3};
p_xi = varCol.nurbs.degree(1);
p_eta = varCol.nurbs.degree(2);
p_zeta = varCol.nurbs.degree(3);
n_zeta = varCol.nurbs.number(3);

element = varCol.element;
weights = varCol.weights;
controlPts = varCol.controlPts;
nurbs = varCol.nurbs;
Upsilon = varCol.Upsilon;

N = varCol.N;
D = varCol.D;
noDofs = varCol.noDofs;

npts = 100;

X = evaluateNURBS(nurbs,[xi,eta,0]);
[R_o, theta, phi] = evaluateProlateCoords(X(1),X(2),X(3),Upsilon);
rArr = linspace(R_o,1.05*R_o,npts);

pArr = zeros(npts,1);
dpArr = zeros(npts,1);

zeta = 0;
% xi = 0.002817541634481;
% eta = 0.004334679437664;
% X = evaluateNURBS(nurbs, [xi,eta,0]);
% [r, theta, phi] = evaluateProlateCoords(X(1),X(2),X(3),Upsilon);
v = zeros(npts,3);
for i = 1:npts
    r = rArr(i);
    X = [sqrt(r^2-Upsilon^2)*sin(theta)*cos(phi);
         sqrt(r^2-Upsilon^2)*sin(theta)*sin(phi);
         r*cos(theta)];
    v(i,:) = X;
    if zeta ~= 1
        parm_pt = pointInversion(nurbs,X,1e-14);
        if abs(parm_pt(3) - 1) < 1e-8
            f = @(parm_pt) norm(X-evaluateNURBS(nurbs, [parm_pt, 1]));
    %         options = optimset('Display','iter','tolx',1e-13);
            options = optimset('tolx',1e-13);
            parm_pt = fminsearchbnd(f,[0.5,0.5], [0,0], [1,1], options);
            zeta = 1;
        else
            zeta = parm_pt(3);  
        end
        xi = parm_pt(1);
        eta = parm_pt(2);
%         zeta = zetaArr(i);
        uniqueXi = unique(Xi);
        uniqueEta = unique(Eta);
        uniqueZeta = unique(Zeta);
        noElementsXi = length(uniqueXi) - 1;
        noElementsEta = length(uniqueEta) - 1;
        noElementsZeta = length(uniqueZeta) - 1;
        xi_idx = findKnotSpan(noElementsXi, 0, xi, uniqueXi);
        eta_idx = findKnotSpan(noElementsEta, 0, eta, uniqueEta);
        zeta_idx = findKnotSpan(noElementsZeta, 0, zeta, uniqueZeta);
        e = xi_idx + noElementsXi*(eta_idx-1) + noElementsXi*noElementsEta*(zeta_idx-1);
        sctr = element(e,:);

        pts = controlPts(sctr,:);
    
        [R, dRdxi, dRdeta, dRdzeta] = NURBS3DBasis(xi, eta, zeta, p_xi, p_eta, p_zeta, Xi, Eta, Zeta, weights);
        J = pts'*[dRdxi' dRdeta' dRdzeta'];
        
        DXDP = dXdP(r,theta,phi,Upsilon);

        e_r = DXDP(:,1)/norm(DXDP(:,1));  

        pArr(i) = R*U(sctr,:);
        dpArr(i) = dot3(J'\[dRdxi; dRdeta; dRdzeta]*U(sctr,:),e_r);
        chi = r;
    else
        R = NURBS3DBasis(xi, eta, 1, p_xi, p_eta, p_zeta, Xi, Eta, Zeta, weights);

        for m = 1:N
            temp = 0;
            temp2 = 0;
            for mt = 1:N
                temp = temp + (1i*k - mt/r)*D(m,mt)*exp(1i*k*(r-chi))*(chi/r)^mt;
                temp2 = temp2 + D(m,mt)*exp(1i*k*(r-chi))*(chi/r)^mt;
            end
            pArr(i) = pArr(i) + temp2*R*U(sctr + noDofs/n_zeta*(m-1),:);
            dpArr(i) = dpArr(i) + temp*R*U(sctr + noDofs/n_zeta*(m-1),:);
        end
    end
    [xi, eta, zeta]
end
if false
    P_0 = varCol.P_0;
    p_exact = P_0*exp(1i*k*(rArr-R_o))*R_o./rArr;
    dp_exact = P_0*(1i*k-1./rArr).*exp(1i*k*(rArr-R_o))*R_o./rArr;
else
    P_0 = varCol.P_0;
    p_exact = scatteredPressureOnRigidSphere2(v(:,1), v(:,2), v(:,3), P_0, k, R_o, 1e-13,varCol.alpha_s_arr,0);
    dp_exact = pressureDeriv_ScatteredPressureOnRigidSphere2(v(:,1), v(:,2), v(:,3), P_0, k, R_o, 1e-13,varCol.alpha_s_arr,0);
    
end

figure(1)
plot(rArr,real(pArr),'*-','color','blue')
hold on
% plot(rArr,real(p_exact),'color','red')
title('p')
hold off
figure(2)
plot(rArr,real(dpArr),'*-','color','blue')
hold on
% plot(rArr,real(dp_exact),'color','red')
title('dp')
hold off

pArr(1)
dpArr(1)
keyboard

X = evaluateNURBS(nurbs,[xi,eta,0]);
normal = [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)];
P_0 = varCol.P_0;
x_0 = [0; 0; 0];
k_vec = -[k*sin(240*pi/180);
          0;
          k*cos(240*pi/180)];
P_inc = P_0*exp(-1i*dot2(k_vec,x_0)).*exp(1i*dot2(k_vec,X'));
dP_inc = 1i*dot2(k_vec,normal).*P_inc;
-dP_inc




