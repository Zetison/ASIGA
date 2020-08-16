function [A, FF, varCol, C, C2] = buildGBEMmatrix(varCol,useSolidDomain)

degree = varCol.degree; % assume p_xi is equal in all patches

index = varCol.index;
noElems = varCol.noElems;
elRange = varCol.elRange;
element = varCol.element;
element2 = varCol.element2;
weights = varCol.weights;
controlPts = varCol.controlPts;
knotVecs = varCol.knotVecs;
pIndex = varCol.pIndex;
patches = varCol.patches;

extraGP = varCol.extraGP;
extraGPBEM = varCol.extraGPBEM;
noDofs = varCol.noDofs;
agpBEM = varCol.agpBEM;
exteriorProblem = varCol.exteriorProblem;
model = varCol.model;
elemsOuter = varCol.elemsOuter;
noDofsInner = varCol.noDofsInner;
noElemsInner = varCol.noElemsInner;
d_p = varCol.patches{1}.nurbs.d_p;

quadMethodBEM = varCol.quadMethodBEM;


Eps = 10*eps;

k = varCol.k;
alpha = 1i/k;
formulation = varCol.formulation;
if strcmp(formulation(end),'C')
    formulation = formulation(1:end-1);
end
switch formulation(2:end)
    case 'BM'
        useCBIE = true;
        useHBIE = true;
        psiType = NaN;
    case 'CBIE'
        useCBIE = true;
        useHBIE = false;
        psiType = NaN;
    case 'HBIE'
        useCBIE = false;
        useHBIE = true;
        psiType = NaN;
    otherwise
        useCBIE = true;
        useHBIE = false;
        psiType = str2double(formulation(end));  
end
useRegul = ~isnan(psiType);

if strcmp(varCol.coreMethod, 'XI')
    useEnrichedBfuns = true;
    d_vec = varCol.d_vec;
else
    useEnrichedBfuns = false;
    d_vec = NaN;
end
if isfield(varCol,'alpha_s')
    no_angles = length(varCol.alpha_s);
else
    no_angles = 1;
end
solveForPtot = varCol.solveForPtot;
if useSolidDomain && ~exteriorProblem
    p_inc = NaN;
    dp_inc = NaN;
    dpdn = @(x,n) 0;
else
    if solveForPtot
        p_inc = varCol.p_inc;
        dp_inc = varCol.dp_inc;
        dpdn = @(x,n) 0;
    else
        p_inc = NaN;
        dp_inc = NaN;
        if strcmp(varCol.BC, 'SHBC')
            dpdn = @(x,n) -varCol.dp_inc(x,n);
        else
            dpdn = varCol.dpdn;
        end
    end
end

if exteriorProblem
    sgn = 1;
else
    sgn = -1;
end

useNeumanProj = varCol.useNeumanProj;
if useNeumanProj
    [U,dU] = projectBC(varCol,SHBC,useCBIE,useHBIE);
else
    U = NaN;
    dU = NaN;
end

[~, ~, diagsMax] = findMaxElementDiameter(patches);
centerPts = findCenterPoints(patches);

n_en = prod(degree+1);

[Q2D_2,W2D_2,Q,W] = getBEMquadData(degree,extraGP,extraGPBEM,quadMethodBEM);
[Q2D,W2D] = gaussTensorQuad(degree+1+extraGP);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_GP = 0;
% nurbs = varCol.nurbs;
% pD = plotBEMGeometry(nurbs,plot_GP,100,1);
% pD = plotBEMGeometry(varCol.nurbs,plot_GP,10,0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

idxRow = zeros(n_en, noElems);
idxRow2 = zeros(n_en, noElems);
idxRow3 = zeros(n_en, noElemsInner);
Avalues = complex(zeros(n_en, noDofs, noElems)); 
d = 3;
if useSolidDomain
    Cvalues = complex(zeros(n_en, d*(noDofs-noDofsInner), noElems)); 
    C2values = complex(zeros(n_en, d*noDofsInner, noElemsInner)); 
else
    Cvalues = NaN; 
    C2values = NaN; 
end
Fvalues = complex(zeros(n_en, noElems, no_angles)); 
totNoQP = 0;
totNoQPnonPolar = 0;

progressBars = varCol.progressBars;
nProgressStepSize = ceil(noElems/1000);
if progressBars
    ppm = ParforProgMon('Building BEM matrix: ', noElems, nProgressStepSize);
else
    ppm = NaN;
end

parfor e_x = 1:noElems
% for e_x = 1:noElems
% for e_x = 1
	if progressBars && mod(e_x,nProgressStepSize) == 0
        ppm.increment();
	end
    if ~plot_GP % to avoid Matlab bug
        pD.plotGP = false;
    end
    patch_x = pIndex(e_x); % New
    knots_x = knotVecs{patch_x}; % New

    Xi_e_x = zeros(d_p,2);
    for i = 1:d_p
        Xi_e_x(i,:) = elRange{i}(index(e_x,i),:);
    end
    J_2_x = prod(Xi_e_x(:,2)-Xi_e_x(:,1))/2^d_p;

    sctr_x = element(e_x,:);
    pts_x = controlPts(sctr_x,:);
    wgts_x = weights(element2(e_x,:)); % New 

    F_e = zeros(n_en, no_angles);
    A_e_temp = zeros(n_en, n_en, noElems);
    if useSolidDomain
        C_e_temp = zeros(n_en, d*n_en, noElems-noElemsInner);
        C2_e_temp = zeros(n_en, d*n_en, noElemsInner);
    end
    
    for gp_x = 1:size(W2D,1)
        pt_x = Q2D(gp_x,:);
        wt_x = W2D(gp_x);

        xi = parent2ParametricSpace(Xi_e_x, pt_x);
        I = findKnotSpans(degree, xi(1,:), knots_x);
        R = NURBSbasis(I, xi, degree, knots_x, wgts_x);
        R_x = R{1};
        dR_xdxi = R{2};
        dR_xdeta = R{3};

        J = [dR_xdxi; dR_xdeta]*pts_x;
        m_1 = J(1,:);
        m_2 = J(2,:);
        crossProd = cross(m_1,m_2);
        J_1_x = norm(crossProd);
        nx = crossProd/J_1_x;

        x = R_x*pts_x;
        fact_x = J_1_x*J_2_x*wt_x;
    
        if useHBIE || useRegul
            h_xi = norm(m_1);
            h_eta = norm(m_2);
            e_xi = m_1/h_xi;
            e_eta = m_2/h_eta;

            v_1 = e_xi;
            v_2 = cross(nx,v_1);
            cosT = dot(e_xi,e_eta);
            sinT = dot(v_2,e_eta);
            dXIdv = [1/h_xi, 0; -cosT/sinT/h_xi, 1/h_eta/sinT];
        else
            dXIdv = NaN;
            v_1 = NaN;
            v_2 = NaN;
        end
        [constants, integrals] = initializeBIE(psiType,useRegul,x,nx,k,model);
        FF_temp = complex(zeros(1, no_angles));
        idxCol = zeros(n_en, noElems);
        idxCol2 = zeros(d*n_en, noElems-noElemsInner);
        idxCol3 = zeros(d*n_en, noElemsInner);
    
        if plot_GP
            pD = plotGP(pD,x,'blue');
        end
        for e_y = 1:noElems   
            [BIE, integrals, FF_temp, sctr_y, noGp, collocationPointIsInElement, pD,~,R_y,r,fact_y,ny] ...
                    = getBEMquadPts(e_y,Q2D_2,W2D_2,Q,W,integrals,FF_temp,...
                    useEnrichedBfuns,k,d_vec,useNeumanProj,solveForPtot,useCBIE,useHBIE,dpdn,U,...
                    x,nx,pt_x(1),pt_x(2),e_x,constants,psiType,useRegul,...
                    degree,pIndex,knotVecs,index,elRange,element,element2,controlPts,weights,...
                    patches,Eps,diagsMax,centerPts,agpBEM,quadMethodBEM,pD);
            idxCol(:,e_y) = sctr_y;
            for i = 1:d
                if ismember(e_y,elemsOuter)
                    idxCol2(i:d:end,e_y-noElemsInner) = d*(sctr_y-1)+i-d*noDofsInner;
                else
                    idxCol3(i:d:end,e_y) = d*(sctr_y-1)+i;
                end
            end
            A_e_temp(:,:,e_y) = A_e_temp(:,:,e_y) + R_x.'*BIE*fact_x;
            if useSolidDomain
                Phi_kTemp = Phi_k(r,k);
                nyR_y = zeros(numel(fact_y),d*n_en);
                for i = 1:d
                    nyR_y(:,i:d:end) = R_y.*ny(:,i);
                end
                if ismember(e_y,elemsOuter)
                    C_e_temp(:,:,e_y-noElemsInner) = C_e_temp(:,:,e_y-noElemsInner) + R_x.'*sum((Phi_kTemp.*fact_y).*nyR_y)*fact_x;
                else
                    C2_e_temp(:,:,e_y) = C2_e_temp(:,:,e_y) + R_x.'*sum((Phi_kTemp.*fact_y).*nyR_y)*fact_x;
                end
            end
            if ~collocationPointIsInElement
                totNoQPnonPolar = totNoQPnonPolar + noGp;
            end
            totNoQP = totNoQP + noGp;
        end
        R_xScaled = getR_x_Coeff(R_x,useEnrichedBfuns,k,d_vec,x,useRegul,integrals,sgn,constants,...
                        psiType,useCBIE,useHBIE,dXIdv,dR_xdxi,dR_xdeta,v_1,v_2,alpha);
             
        A_e_temp(:,:,e_x) = A_e_temp(:,:,e_x) + R_x.'*R_xScaled*fact_x;
        if ~(useSolidDomain && ~exteriorProblem)
            F_eTemp = getF_eTemp(FF_temp,useNeumanProj,solveForPtot,psiType,useCBIE,useHBIE,useRegul,R_x,sctr_x,x,nx,...
                        U,dU,p_inc,dp_inc,dpdn,alpha,integrals,k,constants,sgn);
            F_e = F_e + R_x.'*F_eTemp*fact_x;
        end
    end
    
    idxRow(:,e_x) = sctr_x.';
    if ismember(e_x,elemsOuter)
        idxRow2(:,e_x) = sctr_x.'-noDofsInner;
    else
        idxRow3(:,e_x) = sctr_x.';
    end
    Avalues(:,:,e_x) = matrixAssembly(A_e_temp, idxCol, n_en, noDofs, noElems, 1);
    if useSolidDomain
        if ismember(e_x,elemsOuter)
            Cvalues(:,:,e_x) = matrixAssembly(C_e_temp, idxCol2, n_en, d*(noDofs-noDofsInner), noElems-noElemsInner, 3);
        else
            C2values(:,:,e_x) = matrixAssembly(C2_e_temp, idxCol3, n_en, d*noDofsInner, noElemsInner, 3);
        end
    else
        if ismember(e_x,elemsOuter)
            Cvalues(:,:,e_x) = NaN;
        else
            C2values(:,:,e_x) = NaN;
        end
    end
    Fvalues(:,e_x,:) = F_e;
end
A = matrixAssembly(Avalues, idxRow, n_en, noDofs, noElems, 2);
if useSolidDomain
    C = matrixAssembly(Cvalues(:,:,noElemsInner+1:end), idxRow2(:,noElemsInner+1:end), n_en, noDofs-noDofsInner, noElems-noElemsInner, 4);
    C2 = matrixAssembly(C2values, idxRow3, n_en, noDofsInner, noElemsInner, 4);
else
    C = NaN;
    C2 = NaN;
end
FF = zeros(noDofs,no_angles);
parfor i = 1:no_angles
    FF(:,i) = vectorAssembly(Fvalues(:,:,i), idxRow, noDofs);
end


varCol.totNoQPnonPolar = totNoQPnonPolar;
varCol.totNoQP = totNoQP;
