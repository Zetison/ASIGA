function [A, FF, varCol] = buildBEMmatrixVec(varCol)

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
noElemsPatch = varCol.noElemsPatch;
noPatches = varCol.noPatches;
dofsToRemove = varCol.dofsToRemove;
noDofs = varCol.noDofs;

extraGP = varCol.extraGP;
extraGPBEM = varCol.extraGPBEM;
agpBEM = varCol.agpBEM;
exteriorProblem = varCol.exteriorProblem;
model = varCol.model;

quadMethodBEM = varCol.quadMethodBEM;

Eps = 1e4*eps;

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

SHBC = strcmp(varCol.BC, 'SHBC');
if SHBC
    no_angles = length(varCol.alpha_s);
    p_inc = varCol.p_inc;
    dp_inc = varCol.dp_inc;
else
    no_angles = 1;
    p_inc = NaN;
    dp_inc = NaN;
end
dpdn = varCol.dpdn;

if exteriorProblem
    sgn = 1;
else
    sgn = -1;
end

%% Create collocation points
colBEM_C0 = varCol.colBEM_C0;
if p_xi == 1 && p_eta == 1
    eps_greville_xi = colBEM_C0/(2*p_xi);
    eps_greville_eta = colBEM_C0/(2*p_eta);
else
    eps_greville_xi = colBEM_C0/p_xi;
    eps_greville_eta = colBEM_C0/p_eta;
end
n_cp = noDofs - length(dofsToRemove);
counter2 = 1;
counter = 1;
cp_p = zeros(n_cp,2);
patchIdx = zeros(n_cp,1);
patches = varCol.patches;

[~, ~, diagsMax] = findMaxElementDiameter(patches);
centerPts = findCenterPoints(patches);

for patch = 1:noPatches
    nurbs = patches{patch}.nurbs;
    n_xi = nurbs.number(1);
    n_eta = nurbs.number(2);
    Xi_y = nurbs.knots{1};
    Eta_y = nurbs.knots{2};
    
    if 1
        for j = 1:n_eta
            eta_bar = sum(Eta_y(j+1:j+p_eta))/p_eta;
            if Eta_y(j+1) == Eta_y(j+p_eta)
                if Eta_y(j+1) == Eta_y(j+p_eta+1)
                    eta_bar = eta_bar - eps_greville_eta*(Eta_y(j+p_eta+1)-Eta_y(j));
                else
                    eta_bar = eta_bar + eps_greville_eta*(Eta_y(j+p_eta+1)-Eta_y(j+1));
                end
            end
            for i = 1:n_xi
                if ~any(dofsToRemove == counter)
                    xi_bar = sum(Xi_y(i+1:i+p_xi))/p_xi;
                    if Xi_y(i+1) == Xi_y(i+p_xi)
                        if Xi_y(i+1) == Xi_y(i+p_xi+1)
                            xi_bar = xi_bar - eps_greville_xi*(Xi_y(i+p_xi+1)-Xi_y(i));
                        else
                            xi_bar = xi_bar + eps_greville_xi*(Xi_y(i+p_xi+1)-Xi_y(i+1));
                        end
                    end

                    cp_p(counter2,:) = [xi_bar, eta_bar];
                    patchIdx(counter2) = patch;
                    counter2 = counter2 + 1;
                end
                counter = counter + 1;
            end
        end
    else
        cg_xi = splinesGL(Xi_y,p_xi);
        cg_eta = splinesGL(Eta_y,p_eta);
%         cg_xi = CauchyGalerkin(p_xi, n_xi, Xi);
%         cg_eta = CauchyGalerkin(p_eta, n_eta, Eta);
        for j = 1:n_eta
            eta_bar = cg_eta(j);
            for i = 1:n_xi
                if ~any(dofsToRemove == counter)
                    xi_bar = cg_xi(i);
                    cp_p(counter2,:) = [xi_bar, eta_bar];
                    patchIdx(counter2) = patch;
                    counter2 = counter2 + 1;
                end
                counter = counter + 1;
            end
        end
    end
end
useNeumanProj = varCol.useNeumanProj;
if useNeumanProj
    [U,dU] = projectBC(varCol,SHBC,useCBIE,useHBIE);
else
    U = NaN;
    dU = NaN;
end
eNeighbour = NaN; % to avoid transparency "bug"
createElementTopology

n_en = (p_xi+1)*(p_eta+1);

[Q2D_2,W2D_2,Q,W] = getBEMquadData(p_xi,p_eta,extraGP,extraGPBEM,quadMethodBEM);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_GP = 0;
% pD = plotBEMGeometry(patches,plot_GP,10,0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% A = complex(zeros(1000, noDofs));
% FF = complex(zeros(1000, no_angles));
A = complex(zeros(n_cp, noDofs));
FF = complex(zeros(n_cp, no_angles));
totNoQP = 0;
totNoQPprev = 0;
% for i = 3:n_cp
% for i = 1:n_cp
% for i = 650:n_cp
parfor i = 1:n_cp
% for i = [354,317,319,392]
% for i = 3207:n_cp % BCA rear sail
% for i = 3456:n_cp %  BCA side sail p = 2
% for i = 14472:n_cp
%     i
% for i = 1:n_cp
%     totArea = 0;
    if ~plot_GP % to avoid Matlab bug
        pD.plotGP = false;
    end
    patch = patchIdx(i);
    Xi_x = knotVecs{patch}{1}; % New
    Eta_x = knotVecs{patch}{2}; % New
    uniqueXi = unique(Xi_x);
    uniqueEta = unique(Eta_x);
    noElementsXi = length(uniqueXi)-1;
    noElementsEta = length(uniqueEta)-1;
    
    A_row = complex(zeros(1, noDofs));
    xi_x = cp_p(i,1);
    eta_x = cp_p(i,2);

    xi_idx = findKnotSpan(noElementsXi, 0, xi_x, uniqueXi);
    eta_idx = findKnotSpan(noElementsEta, 0, eta_x, uniqueEta);
    e_x = sum(noElemsPatch(1:patch-1)) + xi_idx + noElementsXi*(eta_idx-1);
    
    idXi_x = index(e_x,1);
    idEta_x = index(e_x,2);

    Xi_e_x = elRangeXi(idXi_x,:);
    Eta_e_x = elRangeEta(idEta_x,:);
    if plot_GP
        if numel(pD.h) > 1
            delete(pD.h(2:end))
        end
        if i == 354
            if strcmp(quadMethodBEM,'Simpson')
                xi_x = parent2ParametricSpace(Xi_e_x, Q(1,1));
                eta_x = parent2ParametricSpace(Eta_e_x, Q(1,1));
            else
                xi_x = parent2ParametricSpace(Xi_e_x, Q{p_xi+1+extraGP}(1,1));
                eta_x = parent2ParametricSpace(Eta_e_x, Q{p_xi+1+extraGP}(1,1));
            end
        end
    end
    
    sctr_x = element(e_x,:);
    pts_x = controlPts(sctr_x,:);
    wgts_x = weights(element2(e_x,:),:); % New

    if useHBIE || useRegul
        singularMapping = true;
        while singularMapping
            [R_x, dR_xdxi, dR_xdeta] = NURBS2DBasis(xi_x, eta_x, p_xi, p_eta, Xi_x, Eta_x, wgts_x);
            x = R_x*pts_x;
            J_temp = [dR_xdxi; dR_xdeta]*pts_x;
            m_1 = J_temp(1,:);
            m_2 = J_temp(2,:);
            crossProd_x = cross(m_1,m_2);
            h_xi = norm(m_1);
            if h_xi < Eps
                eps_greville_eta = 1/(2*p_eta)*(Eta_e_x(2)-Eta_e_x(1));
                if eta_x+eps_greville_eta > Eta_e_x(2)
                    eta_x = eta_x - eps_greville_eta*(Eta_e_x(2)-Eta_e_x(1));
                else
                    eta_x = eta_x + eps_greville_eta*(Eta_e_x(2)-Eta_e_x(1));
                end
                continue
            end
            h_eta = norm(m_2);
            if h_eta < Eps
                eps_greville_xi = 1/(2*p_xi)*(Xi_e_x(2)-Xi_e_x(1));
                if xi_x+eps_greville_xi > Xi_e_x(2)
                    xi_x = xi_x - eps_greville_xi*(Xi_e_x(2)-Xi_e_x(1));
                else
                    xi_x = xi_x + eps_greville_xi*(Xi_e_x(2)-Xi_e_x(1));
                end
                continue
            end
            singularMapping = false;
        end
        e_xi = m_1/h_xi;
        e_eta = m_2/h_eta;

        v_1 = m_1/norm(m_1);
        nx = crossProd_x/norm(crossProd_x);
        v_2 = cross(nx,v_1);
        cosT = dot(e_xi,e_eta);
        sinT = dot(v_2,e_eta);
        dXIdv = [1/h_xi, 0; -cosT/sinT/h_xi, 1/h_eta/sinT];
    else
        R_x = NURBS2DBasis(xi_x, eta_x, p_xi, p_eta, Xi_x, Eta_x, wgts_x);
        dR_xdxi = NaN;
        dR_xdeta = NaN;
        x = R_x*pts_x;
        nx = NaN;       
        dXIdv = NaN;  
        v_1 = NaN;
        v_2 = NaN;
    end
%     if x(1) < -15.9 && x(1) > -17 && x(2) < 1.2 && x(2) > 1.19 && abs(x(3) - 4) < 1e6*eps
%         keyboard
%     else
%         continue
%     end
    [constants, integrals] = initializeBIE(psiType,useRegul,x,nx,k,model);
    
    FF_temp = complex(zeros(1, no_angles));
    if plot_GP
        pD = plotGP(pD,x,'blue');
    end
    [adjacentElements, xi_x_tArr,eta_x_tArr] = getAdjacentElements(e_x,xi_x,eta_x,Xi_e_x,Eta_e_x,eNeighbour,Eps);
    for e_y = 1:noElems  
%         if e_y == 32
%             keyboard
%         end
        [BIE, integrals, FF_temp, sctr_y, noGp, pD] = getBEMquadPts(e_y,Q2D_2,W2D_2,Q,W,integrals,FF_temp,...
                useEnrichedBfuns,k,d_vec,useNeumanProj,SHBC,useCBIE,useHBIE,dpdn,U,...
                x,nx,xi_x_tArr,eta_x_tArr,xi_x,eta_x,adjacentElements,constants,psiType,useRegul,...
                p_xi, p_eta,pIndex,knotVecs,index,elRangeXi,elRangeEta,element,element2,controlPts,weights,...
                patches,Eps,diagsMax,centerPts,agpBEM,quadMethodBEM,pD);
        for j = 1:n_en
            A_row(sctr_y(j)) = A_row(sctr_y(j)) + BIE(j);
        end
        totNoQP = totNoQP + noGp;
    end
    if plot_GP
        figureFullScreen(gcf)
%         totNoQP-totNoQPprev
%         totNoQPprev = totNoQP;
%         export_fig(['../../graphics/BEM/S1_' num2str(i) '_extraGPBEM' num2str(extraGPBEM) '_agpBEM' num2str(agpBEM) '_' quadMethodBEM], '-png', '-transparent', '-r300')
        keyboard
    end
    R_xScaled = getR_x_Coeff(R_x,useEnrichedBfuns,k,d_vec,x,useRegul,integrals,sgn,constants,...
                    psiType,useCBIE,useHBIE,dXIdv,dR_xdxi,dR_xdeta,v_1,v_2,alpha);
    for j = 1:n_en
        A_row(sctr_x(j)) = A_row(sctr_x(j)) + R_xScaled(j);
    end     
    if any(isinf(A_row)) || any(isnan(A_row))
        error(['Problems at i = ' num2str(i)])
    end
    A(i,:) = A_row;
    FF(i,:) = getF_eTemp(FF_temp,useNeumanProj,SHBC,psiType,useCBIE,useHBIE,useRegul,R_x,sctr_x,x,nx,...
                U,dU,p_inc,dp_inc,dpdn,alpha,integrals,k,constants,sgn);
end

% totNoQP
varCol.totNoQP = totNoQP;
