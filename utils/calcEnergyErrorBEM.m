function relError = calcEnergyErrorBEM(varCol)

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

extraGP = varCol.extraGP;
extraGPBEM = varCol.extraGPBEM;
agpBEM = varCol.agpBEM;
analytic = varCol.analytic;
solveForPtot = varCol.solveForPtot;

quadMethodBEM = varCol.quadMethodBEM;
U = varCol{1}.U;

Eps = 10*eps;

kPhi = 1i;
Phi_ = @(r) exp(1i*kPhi*r)./(4*pi*r);

k = varCol.k;

if strcmp(varCol.coreMethod, 'XI')
    useEnrichedBfuns = true;
    d_vec = varCol.d_vec;
else
    useEnrichedBfuns = false;
    d_vec = NaN;
end
[~, ~, diagsMax] = findMaxElementDiameter(patches);
centerPts = findCenterPoints(patches);


[Q2D_2,W2D_2,Q,W] = getBEMquadData(p_xi,p_eta,extraGP,extraGPBEM,quadMethodBEM);
[Q2D,W2D] = tensorQuad(p_xi+1+extraGP,p_eta+1+extraGP);
% [W2D,Q2D] = gaussianQuadNURBS(p_xi+10,p_eta+10);  

dpdn = @(x,n) 0;
if solveForPtot
    analytic = @(x) analytic(x) + varCol.p_inc(x);
end
Error = 0;
normalization = 0;
parfor e_x = 1:noElems
% for e_x = 1:noElems
    patch_x = pIndex(e_x); % New
    Xi_x = knotVecs{patch_x}{1}; % New
    Eta_x = knotVecs{patch_x}{2}; % New

    idXi_x = index(e_x,1);
    idEta_x = index(e_x,2);

    Xi_e_x = elRangeXi(idXi_x,:);
    Eta_e_x = elRangeEta(idEta_x,:);

    sctr_x = element(e_x,:);
    pts_x = controlPts(sctr_x,:);
    wgts_x = weights(element2(e_x,:)); % New 
    
    J_2_x = 0.25*(Xi_e_x(2)-Xi_e_x(1))*(Eta_e_x(2)-Eta_e_x(1));
    
    for gp_x = 1:size(W2D,1)
        pt_x = Q2D(gp_x,:);
        wt_x = W2D(gp_x);

        xi_x  = parent2ParametricSpace(Xi_e_x, pt_x(1));
        eta_x = parent2ParametricSpace(Eta_e_x,pt_x(2));
        [R_x, dR_xdxi, dR_xdeta] = NURBS2DBasis(xi_x, eta_x, p_xi, p_eta, Xi_x, Eta_x, wgts_x);

        J = [dR_xdxi; dR_xdeta]*pts_x;
        m_1 = J(1,:);
        m_2 = J(2,:);
        crossProd = cross(m_1,m_2);
        J_1_x = norm(crossProd);
        nx = crossProd/J_1_x;

        x = R_x*pts_x;
        fact_x = J_1_x*J_2_x*wt_x;
        analytic_x = analytic(x);
        u_x = R_x*U(sctr_x);
        for e_y = 1:noElems   
            [~, ~, ~, sctr_y, ~, ~, ~,y,R_y,r,fact_y] = getBEMquadPts(e_y,Q2D_2,W2D_2,Q,W,NaN,NaN,...
                    useEnrichedBfuns,k,d_vec,0,solveForPtot,1,0,dpdn,NaN,...
                    x,nx,pt_x(1),pt_x(2),e_x,NaN,NaN,0,...
                    p_xi, p_eta,pIndex,knotVecs,index,elRangeXi,elRangeEta,element,element2,controlPts,weights,...
                    patches,Eps,diagsMax,centerPts,agpBEM,quadMethodBEM);
            analytic_y = analytic(y);
            Error = Error + conj(analytic_x-u_x)*(Phi_(r).*(analytic_y-R_y*U(sctr_y))).'*fact_y*fact_x;
            normalization = normalization + conj(analytic_x)*(Phi_(r).*analytic_y).'*fact_y*fact_x;
        end
    end
end

relError = 100*sqrt(abs(Error/normalization));


