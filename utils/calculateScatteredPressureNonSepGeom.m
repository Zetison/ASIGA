function p_h = calculateScatteredPressureNonSepGeom(varCol, U, P_far, computeFarField)
error('Depricated, use calculateScattere instead')

p_xi = varCol.degree(1); % assume p_xi is equal in all patches
p_eta = varCol.degree(2); % assume p_eta is equal in all patches
N = varCol.N;
solveForPtot = varCol.solveForPtot;
SHBC = strcmp(varCol.BC, 'SHBC');

index = varCol.index;
noElems = varCol.noElems;
elRangeXi = varCol.elRange{1};
elRangeEta = varCol.elRange{2};
element = varCol.element;
element2 = varCol.element2;
weights = varCol.weights;
controlPts = varCol.controlPts;
knotVecs = varCol.knotVecs;
pIndex = varCol.pIndex;

dp_inc = varCol.dp_inc;
homNeumanCond = false;
dpdn = @(x,n) -dp_inc(x,n);

Phi_k = varCol.Phi_k;
dPhi_kdny = varCol.dPhi_kdny;

extraGP = varCol.extraGP;

x_0 = varCol.x_0;
A_2 = varCol.A_2;

k = varCol.k;
Upsilon = varCol.Upsilon;

noDofs = varCol.noDofs;

D = varCol.D;
if 0
    X = evaluateNURBS(nurbs,[0,0]);
    chi = evaluateProlateCoords(X(1),X(2),X(3),Upsilon);
    p_h = zeros(1,size(P_far,2), size(P_far,3));
    uniqueXi = unique(Xi);
    uniqueEta = unique(Eta);
    noElementsXi = length(uniqueXi) - 1;
    noElementsEta = length(uniqueEta) - 1;
    for i = 1:size(P_far,2)
        for j = 1:size(P_far,3)
            x_far = P_far(1,i,j);
            y_far = P_far(2,i,j);
            z_far = P_far(3,i,j);
            r_far = evaluateProlateCoords(x_far,y_far,z_far,Upsilon);
            
            [xi, eta] = prolateCoordsToParamCoords(nurbs, x_far, y_far, z_far, Upsilon, chi);

            R = NURBS2DBasis(xi, eta, p_xi, p_eta, Xi, Eta, weights);


            xi_idx = findKnotSpan(noElementsXi, 0, xi, uniqueXi);
            eta_idx = findKnotSpan(noElementsEta, 0, eta, uniqueEta);
            e = xi_idx + noElementsXi*(eta_idx-1);
            sctr = element(e,:);
            
%             vSurf = R*controlPts(sctr,:);
% 
%             x = vSurf(1);
%             y = vSurf(2);
%             z = vSurf(3);
% 
%             chi = evaluateProlateCoords(x,y,z,Upsilon);

            u = 0;
            for m = 1:N
                Q_m = 0;
                for mt = 1:N
                    Q_m = Q_m + (chi/r_far)^mt*D(m,mt);
                end
                phi_m = exp(1i*k*(r_far-chi))*Q_m;

                u = u + phi_m*(R*U(sctr + noDofs*(m-1)));    
            end
        end
        p_h(1,i,:) = u;
    end
        
else
    [W2D,Q2D] = gaussianQuadNURBS(p_xi+1+extraGP,p_eta+1+extraGP);  

    p_h = zeros(size(P_far,1),1);
    
%     for e = 1:noElems %
    parfor e = 1:noElems
        patch = pIndex(e); % New
        Xi = knotVecs{patch}{1}; % New
        Eta = knotVecs{patch}{2}; % New

        idXi = index(e,1);
        idEta = index(e,2);

        Xi_e = elRangeXi(idXi,:);
        Eta_e = elRangeEta(idEta,:);

        sctr = element(e,:);
        pts = controlPts(sctr,:);
        wgts = weights(element2(e,:)); % New   
        
        J_2 = 0.25*(Xi_e(2)-Xi_e(1))*(Eta_e(2)-Eta_e(1));

        for gp = 1:size(W2D,1)
            pt = Q2D(gp,:);
            wt = W2D(gp);

            xi   = parent2ParametricSpace(Xi_e,  pt(1));
            eta  = parent2ParametricSpace(Eta_e, pt(2));

            [R, dRdxi, dRdeta] = NURBS2DBasis(xi, eta, p_xi, p_eta, Xi, Eta, wgts);

            J = pts'*[dRdxi' dRdeta'];
            crossProd = cross(J(:,1),J(:,2));
            J_1 = norm(crossProd);
            n = crossProd.'/J_1;

            Y = R*pts;
            Yt = A_2*(Y-x_0)';
            xt = Yt(1);
            yt = Yt(2);
            zt = Yt(3);

            r_a = evaluateProlateCoords(xt,yt,zt,Upsilon);
%             DPDX = dPdX(xt,yt,zt,Upsilon,r_a,c1,c2)*A_2;

%             J3 = DPDX(2:3,:)*J;

%             [dchidtheta, dchidphi] = dchidP(DPDX(1,:),J,J3);

%             DXDP = A_2\dXdP(r,theta,phi,Upsilon);

%             h_r = norm(DXDP(:,1));
%             h_theta = norm(DXDP(:,2));
%             h_phi = norm(DXDP(:,3));
% 
%             e_r = DXDP(:,1)/h_r;
%             e_theta = DXDP(:,2)/h_theta;
%             e_phi = DXDP(:,3)/h_phi;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %         p_h = R_fun*U(sctr,:,:);
    %         dp_h = dot3(matrix3Dprod(J'\[dRdxi; dRdeta; dRdzeta],U(sctr,:,:)),normal);
            p_h_gp = zeros(1,size(U,2));
            for m = 1:N
                temp = 0;
                temp2 = 0;
                for mt = 1:N
                    temp = temp + (1i*k - mt/r_a)*D(m,mt);
                    temp2 = temp2 + D(m,mt);
                end
                p_h_gp = p_h_gp + temp2*R*U(sctr + noDofs*(m-1),:);
            end    
            xmy = -elementAddition(Y, -P_far);

            r = norm2(xmy);
            
            if computeFarField
                x_d_n = dot3(P_far, n.')./norm2(P_far);
                x_d_y = dot3(P_far, Y.')./norm2(P_far);
                p_h = p_h + -1/(4*pi)*1i*k* (p_h_gp.').*x_d_n.*exp(-1i*k*x_d_y)* J_1 * J_2 * wt;  
            else
                p_h = p_h + (p_h_gp.').*dPhi_kdny(xmy,r,n,k)* J_1 * J_2 * wt;  
            end
            if ~homNeumanCond
                dp_h_gp = dpdn(Y,n).';
                if computeFarField
                    p_h = p_h + -1/(4*pi)*dp_h_gp.*exp(-1i*k*x_d_y)* J_1 * J_2 * wt;  
                else
                    p_h = p_h + -dp_h_gp.*Phi_k(r,k)* J_1 * J_2 * wt;  
                end
            end
        end
    end
end

