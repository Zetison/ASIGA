function [L2Error,H1error,H1sError] = calcError2D_FEM(varCol, U, exactNormalizationH1s)

%% Extract all needed data from options and varCol
p_xi = varCol.patches{1}.nurbs.degree(1);
p_eta = varCol.patches{1}.nurbs.degree(2);
controlPts = varCol.patches{1}.controlPts;

j_lag_xi = createIdxMap(p_xi+1);
i_lag_xi = repmat((1:p_xi+1)',1,p_xi);
j_lag_eta = createIdxMap(p_eta+1);
i_lag_eta = repmat((1:p_eta+1)',1,p_eta);
% GLL_p = GLLpoints(p+1);  
% GLL_q = GLLpoints(q+1);  
GLL_xi = linspace(-1,1,p_xi+1);  
GLL_eta = linspace(-1,1,p_eta+1);  

noElems = varCol.patches{1}.noElems;
element = varCol.patches{1}.element;

analytic = varCol.analytic;
danalytic = varCol.danalytic;

%% Preallocation and initiallizations
[W2D,Q2D] = gaussianQuadNURBS(p_xi+5,p_eta+5); 

%% Build global matrices
% warning('Not running in Parallel')
% keyboard
L2Error = 0;
H1error = 0;
H1sError = 0;
normalizationL2 = 0;
normalizationH1 = 0;
normalizationH1s = 0;
parfor e = 1:noElems
% for e = 1:noElems    
    sctr = element(e,:);
    pts = controlPts(sctr,:);
    U_sctr = U(sctr);
    
    for gp = 1:size(W2D,1)
        pt = Q2D(gp,:);
        wt = W2D(gp);
        
        [L_xi,dL_xi] = lagrangeBasis(pt(1),i_lag_xi,GLL_xi,j_lag_xi);
        [L_eta,dL_eta] = lagrangeBasis(pt(2),i_lag_eta,GLL_eta,j_lag_eta);
        
        R = kron(L_xi,L_eta);
        dRdxi = kron(dL_xi,L_eta);
        dRdeta = kron(L_xi,dL_eta);
        J = pts'*[dRdxi' dRdeta'];
        J_1 = det(J);
        dRdX = J'\[dRdxi; dRdeta];
        
        u_h = R*U_sctr;
        du_h = dRdX*U_sctr;
        y = R*pts;
        u = analytic(y(1),y(2));
        du = danalytic(y(1),y(2));
        L2Error = L2Error + abs(u-u_h)^2 * abs(J_1) * wt;
        H1error = H1error + (norm(du-du_h)^2 + abs(u-u_h)^2) * abs(J_1) * wt;
        H1sError = H1sError + norm(du-du_h)^2 * abs(J_1) * wt;
        normalizationL2 = normalizationL2 + abs(u)^2 * abs(J_1) * wt;
        normalizationH1 = normalizationH1 + (norm(du)^2 + abs(u)^2) * abs(J_1) * wt;
        normalizationH1s = normalizationH1s + norm(du)^2 * abs(J_1) * wt;
    end
end

L2Error = 100*sqrt(L2Error/normalizationL2);
H1error = 100*sqrt(H1error/normalizationH1);
if nargin > 2 
    H1sError = 100*sqrt(H1sError)/exactNormalizationH1s;
else
    H1sError = 100*sqrt(H1sError/normalizationH1s);
end
