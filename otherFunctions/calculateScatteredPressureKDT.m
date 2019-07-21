function p_h = calculateScatteredPressureKDT(varCol, P_far, plotFarField)


if ~strcmp(varCol.BC, 'SHBC')
    error('This is not implemented')
end
p_xi = varCol.degree(1); % assume p_xi is equal in all patches
p_eta = varCol.degree(2); % assume p_eta is equal in all patches

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
patches = varCol.patches;
agpBEM = varCol.agpBEM;
extraGPBEM = varCol.extraGPBEM;
formulation = varCol.formulation;
Eps = 10*eps;

d_vec = varCol.d_vec;
k = varCol.k;
[~, ~, diagsMax] = findMaxElementDiameter(patches);
centerPts = findCenterPoints(patches);

p_inc = varCol.p_inc;

extraGP = varCol.extraGP;
[~,W2D_2,Q,W] = getBEMquadData(p_xi,p_eta,extraGP,extraGPBEM,'adaptive');

[W2D,Q2D] = gaussianQuadNURBS(p_xi+1+extraGP,p_eta+1+extraGP);  
% [W2D,Q2D] = gaussianQuadNURBS(20,20);
p_h = zeros(size(P_far,1),length(k));
switch formulation
    case 'MS1'
        for e_y = 1:noElems
%         parfor e = 1:noElems
            patch_y = pIndex(e_y); % New
            Xi_y = knotVecs{patch_y}{1}; % New
            Eta_y = knotVecs{patch_y}{2}; % New

            idXi_y = index(e_y,1);
            idEta_y = index(e_y,2);

            Xi_e_y = elRangeXi(idXi_y,:);
            Eta_e_y = elRangeEta(idEta_y,:);

            sctr_y = element(e_y,:);
            pts_y = controlPts(sctr_y,:);
            wgts_y = weights(element2(e_y,:)); % New   

            J_2_y = 0.25*(Xi_e_y(2)-Xi_e_y(1))*(Eta_e_y(2)-Eta_e_y(1));

            for gp_y = 1:size(W2D,1)
                pt_y = Q2D(gp_y,:);
                wt_y = W2D(gp_y);

                xi_y   = parent2ParametricSpace(Xi_e_y,  pt_y(1));
                eta_y  = parent2ParametricSpace(Eta_e_y, pt_y(2));

                [R_y, dRdxi_y, dRdeta_y] = NURBS2DBasis(xi_y, eta_y, p_xi, p_eta, Xi_y, Eta_y, wgts_y);

                J_y = pts_y'*[dRdxi_y' dRdeta_y'];
                crossProd_y = cross(J_y(:,1),J_y(:,2));
                J_1_y = norm(crossProd_y);
                n_y = crossProd_y/J_1_y;

                y = R_y*pts_y;
                fact_y = J_1_y * J_2_y * wt_y;
                y_d_n = P_far*n_y./norm2(P_far);
                x_d_y = P_far*y.'./norm2(P_far);
                dPhiFarField = -1i*k*y_d_n.*exp(-1i*k*x_d_y)/(4*pi);
                
                for e_z = 1:noElems   
                    e_z
                    [z,n_z,fact_z] = getBEMquadPtsKDT(e_y,e_z,W2D_2,Q,W,y,pt_y(1),pt_y(2),...
                        p_xi, p_eta,pIndex,knotVecs,index,elRangeXi,elRangeEta,element,element2,controlPts,weights,...
                        patches,Eps,diagsMax,centerPts,agpBEM);
                    p_h_gp = 4*p_inc(z);
                    p_h_gp(or(n_z*d_vec > 0,repmat(sum(n_z.*(y-z) > 0,2),1,size(d_vec,2)))) = 0;
                    ymz = repmat(y,size(z,1),1) - z;
                    r = norm2(ymz);

                    if plotFarField
                        p_h = p_h + (dPhiFarField * (dPhi_kdny(ymz,r,n_y.',k).* fact_z).'*p_h_gp).' * fact_y;  
                    else
%                         xmy = P_far - repmat(x,size(P_far,1),1);
%                         r = norm2(xmy);
%                         p_h = p_h + p_h_gp.'.*dPhi_kdny(xmy,r,n_x.',k)* fact_x;  
                    end
                end
            end
        end
        
    otherwise
        % for e = 1:noElems
        parfor e_y = 1:noElems
            patch_y = pIndex(e_y); % New
            Xi_y = knotVecs{patch_y}{1}; % New
            Eta_y = knotVecs{patch_y}{2}; % New

            idXi_y = index(e_y,1);
            idEta_y = index(e_y,2);

            Xi_e_y = elRangeXi(idXi_y,:);
            Eta_e_y = elRangeEta(idEta_y,:);

            sctr_y = element(e_y,:);
            pts_y = controlPts(sctr_y,:);
            wgts_y = weights(element2(e_y,:)); % New   

            J_2_y = 0.25*(Xi_e_y(2)-Xi_e_y(1))*(Eta_e_y(2)-Eta_e_y(1));

            for gp_y = 1:size(W2D,1)
                pt_y = Q2D(gp_y,:);
                wt_y = W2D(gp_y);

                xi_y   = parent2ParametricSpace(Xi_e_y,  pt_y(1));
                eta_y  = parent2ParametricSpace(Eta_e_y, pt_y(2));

                [R_y, dRdxi_y, dRdeta_y] = NURBS2DBasis(xi_y, eta_y, p_xi, p_eta, Xi_y, Eta_y, wgts_y);

                J_y = pts_y'*[dRdxi_y' dRdeta_y'];
                crossProd_y = cross(J_y(:,1),J_y(:,2));
                J_1_y = norm(crossProd_y);
                n_y = crossProd_y/J_1_y;

                y = R_y*pts_y;

                p_h_gp = 2*p_inc(y);
                p_h_gp(n_y.'*d_vec > 0) = 0;

                if plotFarField
                    y_d_n = P_far*n_y./norm2(P_far);
                    x_d_y = P_far*y.'./norm2(P_far);
                    p_h = p_h - 1/(4*pi)*1i*k.*p_h_gp.'.*y_d_n.*exp(-1i*k*x_d_y)* J_1_y * J_2_y * wt_y;  
                else
                    xmy = P_far - repmat(y,size(P_far,1),1);
                    r = norm2(xmy);
                    p_h = p_h + p_h_gp.'.*dPhi_kdny(xmy,r,n_y.',k)* J_1_y * J_2_y * wt_y;  
                end
            end
        end
end


