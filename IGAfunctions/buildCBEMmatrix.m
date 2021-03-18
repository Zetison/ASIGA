function task = buildCBEMmatrix(task)

degree = task.varCol{1}.degree; % assume degree is equal in all patches

index = task.varCol{1}.index;
noElems = task.varCol{1}.noElems;
elRange = task.varCol{1}.elRange;
element = task.varCol{1}.element;
element2 = task.varCol{1}.element2;
weights = task.varCol{1}.weights;
controlPts = task.varCol{1}.controlPts;
knotVecs = task.varCol{1}.knotVecs;
pIndex = task.varCol{1}.pIndex;
noElemsPatch = task.varCol{1}.noElemsPatch;
noPatches = task.varCol{1}.noPatches;
dofsToRemove = task.varCol{1}.dofsToRemove;
noDofs = task.varCol{1}.noDofs;

extraGP = task.misc.extraGP;
extraGPBEM = task.bem.extraGPBEM;
agpBEM = task.bem.agpBEM;
exteriorProblem = task.misc.exteriorProblem;
model = task.misc.model;
colMethod = task.bem.colMethod;

quadMethodBEM = task.bem.quadMethodBEM;

Eps = 1e4*eps;

k = task.varCol{1}.k;
alpha = 1i/k;

formulation = task.misc.formulation;
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

if strcmp(task.misc.coreMethod, 'XI')
    useEnrichedBfuns = true;
    d_vec = task.varCol{1}.d_vec;
else
    useEnrichedBfuns = false;
    d_vec = NaN;
end

SHBC = strcmp(task.misc.BC, 'SHBC');
if SHBC
    no_angles = length(task.ffp.alpha_s);
else
    no_angles = 1;
end
solveForPtot = task.misc.solveForPtot;
if solveForPtot
    p_inc = task.p_inc_;
    dp_inc = task.dp_inc_;
    dpdn = @(x,n) 0;
else
    p_inc = NaN;
    dp_inc = task.dp_inc_;
    if SHBC
        dpdn = @(x,n) -dp_inc_(x,n);
    else
        dpdn = task.dpdn_;
    end
end

if exteriorProblem
    sgn = 1;
else
    sgn = -1;
end

%% Create collocation points
colBEM_C0 = task.bem.colBEM_C0;
if all(degree == 1)
    eps_greville = colBEM_C0./(2*degree);
else
    eps_greville = colBEM_C0./degree;
end
n_cp = noDofs - length(dofsToRemove);
counter2 = 1;
counter = 1;
cp_p = zeros(n_cp,2);
patchIdx = zeros(n_cp,1);
patches = task.varCol{1}.patches;

[~, ~, diagsMax] = findMaxElementDiameter(patches);
centerPts = findCenterPoints(patches);

for patch = 1:noPatches
    nurbs = patches{patch}.nurbs;
    n_xi = nurbs.number(1);
    n_eta = nurbs.number(2);
    Xi_y = nurbs.knots{1};
    Eta_y = nurbs.knots{2};
    
    switch colMethod
        case 'Grev'
            for j = 1:n_eta
                eta_bar = sum(Eta_y(j+1:j+degree(2)))/degree(2);
                if Eta_y(j+1) == Eta_y(j+degree(2))
                    if Eta_y(j+1) == Eta_y(j+degree(2)+1)
                        eta_bar = eta_bar - eps_greville(2)*(Eta_y(j+degree(2)+1)-Eta_y(j));
                    else
                        eta_bar = eta_bar + eps_greville(2)*(Eta_y(j+degree(2)+1)-Eta_y(j+1));
                    end
                end
                for i = 1:n_xi
                    if ~any(dofsToRemove == counter)
                        xi_bar = sum(Xi_y(i+1:i+degree(1)))/degree(1);
                        if Xi_y(i+1) == Xi_y(i+degree(1))
                            if Xi_y(i+1) == Xi_y(i+degree(1)+1)
                                xi_bar = xi_bar - eps_greville(1)*(Xi_y(i+degree(1)+1)-Xi_y(i));
                            else
                                xi_bar = xi_bar + eps_greville(1)*(Xi_y(i+degree(1)+1)-Xi_y(i+1));
                            end
                        end

                        cp_p(counter2,:) = [xi_bar, eta_bar];
                        patchIdx(counter2) = patch;
                        counter2 = counter2 + 1;
                    end
                    counter = counter + 1;
                end
            end
        case 'CG'
            cg_xi = CauchyGalerkin(degree(1), n_xi, Xi_y);
            cg_eta = CauchyGalerkin(degree(2), n_eta, Eta_y);
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
        case 'GL'
            cg_xi = splinesGL(Xi_y,degree(1));
            cg_eta = splinesGL(Eta_y,degree(2));
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
useNeumanProj = task.bem.useNeumanProj;
if useNeumanProj
    [U,dU] = projectBC(task.varCol{1},SHBC,useCBIE,useHBIE);
else
    U = NaN;
    dU = NaN;
end
eNeighbour = createElementTopology(task);

n_en = prod(degree+1);

[Q2D_2,W2D_2,Q,W] = getBEMquadData(degree,extraGP(1:2),extraGPBEM,quadMethodBEM);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_GP = 0;
% nurbs = task.varCol{1}.nurbs;
% pD = plotBEMGeometry(nurbs,plot_GP,100,1);
% pD = plotBEMGeometry(task.varCol{1}.nurbs,plot_GP,10,0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% A = complex(zeros(1000, noDofs));
% FF = complex(zeros(1000, no_angles));
A = complex(zeros(n_cp, noDofs));
FF = complex(zeros(n_cp, no_angles));
totNoQPnonPolar = 0;
totNoQP = 0;
totNoQPprev = 0;

nProgressStepSize = ceil(n_cp/1000);
progressBars = task.misc.progressBars;
if progressBars
    ppm = ParforProgMon('Building BEM matrix: ', n_cp, nProgressStepSize);
else
    ppm = NaN;
end

% for i = 247
% for i = 1:n_cp
parfor i = 1:n_cp
	if progressBars && mod(i,nProgressStepSize) == 0
        ppm.increment();
	end
%     totArea = 0;
    if ~plot_GP % to avoid Matlab bug
        pD.plotGP = false;
    end
    patch = patchIdx(i);
    knots = knotVecs{patch};
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

    Xi_e_x = elRange{1}(idXi_x,:);
    Eta_e_x = elRange{2}(idEta_x,:);
    if plot_GP
        if numel(pD.h) > 1
            delete(pD.h(2:end))
        end
        if i == 354
            if strcmp(quadMethodBEM,'Simpson')
                xi_x = parent2ParametricSpace(Xi_e_x, Q(1,1));
                eta_x = parent2ParametricSpace(Eta_e_x, Q(1,1));
            else
                xi_x = parent2ParametricSpace(Xi_e_x, Q{degree(1)+1+extraGP(1)}(1,1));
                eta_x = parent2ParametricSpace(Eta_e_x, Q{degree(1)+1+extraGP(1)}(1,1));
            end
        end
    end
    
    sctr_x = element(e_x,:);
    pts_x = controlPts(sctr_x,:);
    wgts_x = weights(element2(e_x,:),:); % New

    if useHBIE || useRegul
        singularMapping = true;
        eps_greville = zeros(1,2);
        while singularMapping
            xi = [xi_x, eta_x];
            I = findKnotSpans(degree, xi(1,:), knots);
            R = NURBSbasis(I, xi, degree, knots, wgts_x);
            R_x = R{1};
            dR_xdxi = R{2};
            dR_xdeta = R{3};
            x = R_x*pts_x;
            J_temp = [dR_xdxi; dR_xdeta]*pts_x;
            m_1 = J_temp(1,:);
            m_2 = J_temp(2,:);
            crossProd_x = cross(m_1,m_2);
            h_xi = norm(m_1);
            if h_xi < Eps
                eps_greville(2) = 1/(2*degree(2))*(Eta_e_x(2)-Eta_e_x(1));
                if eta_x+eps_greville(2) > Eta_e_x(2)
                    eta_x = eta_x - eps_greville(2)*(Eta_e_x(2)-Eta_e_x(1));
                else
                    eta_x = eta_x + eps_greville(2)*(Eta_e_x(2)-Eta_e_x(1));
                end
                continue
            end
            h_eta = norm(m_2);
            if h_eta < Eps
                eps_greville(1) = 1/(2*degree(1))*(Xi_e_x(2)-Xi_e_x(1));
                if xi_x+eps_greville(1) > Xi_e_x(2)
                    xi_x = xi_x - eps_greville(1)*(Xi_e_x(2)-Xi_e_x(1));
                else
                    xi_x = xi_x + eps_greville(1)*(Xi_e_x(2)-Xi_e_x(1));
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
        xi = [xi_x, eta_x];
        I = findKnotSpans(degree, xi(1,:), knots);
        R = NURBSbasis(I, xi, degree, knots, wgts_x);
        R_x = R{1};
        dR_xdxi = NaN;
        dR_xdeta = NaN;
        x = R_x*pts_x;
        nx = NaN;       
        dXIdv = NaN;  
        v_1 = NaN;
        v_2 = NaN;
    end
    [constants, integrals] = initializeBIE(psiType,useRegul,x,nx,k,model);
    
    FF_temp = complex(zeros(1, no_angles));
    if plot_GP
        pD = plotGP(pD,x,'blue');
    end
    [adjacentElements, xi_x_tArr,eta_x_tArr] = getAdjacentElements(e_x,xi_x,eta_x,Xi_e_x,Eta_e_x,eNeighbour,Eps);
    for e_y = 1:noElems  
        [BIE, integrals, FF_temp, sctr_y, noGp, collocationPointIsInElement, pD] = getBEMquadPts(e_y,Q2D_2,W2D_2,Q,W,integrals,FF_temp,...
                useEnrichedBfuns,k,d_vec,useNeumanProj,solveForPtot,useCBIE,useHBIE,dpdn,U,...
                x,nx,xi_x_tArr,eta_x_tArr,adjacentElements,constants,psiType,useRegul,...
                degree,pIndex,knotVecs,index,elRange,element,element2,controlPts,weights,...
                patches,Eps,diagsMax,centerPts,agpBEM,quadMethodBEM,pD);
        for j = 1:n_en
            A_row(sctr_y(j)) = A_row(sctr_y(j)) + BIE(j);
        end
        if ~collocationPointIsInElement
            totNoQPnonPolar = totNoQPnonPolar + noGp;
        end
        totNoQP = totNoQP + noGp;
    end
%     rms(A_row)
    if plot_GP
%         figureFullScreen(gcf)
%         totNoQP-totNoQPprev
%         totNoQPprev = totNoQP;
%         export_fig(['../../graphics/BEM/S1_' num2str(i) '_extraGPBEM' num2str(extraGPBEM) '_agpBEM' num2str(agpBEM) '_' quadMethodBEM], '-png', '-transparent', '-r300')
%         export_fig(['../../graphics/BEM/S1_' num2str(i) '_extraGPBEM' num2str(extraGPBEM) '_agpBEM' num2str(agpBEM) '_' quadMethodBEM], '-png', '-transparent', '-r200')
%         keyboard
    end
    R_xScaled = getR_x_Coeff(R_x,useEnrichedBfuns,k,d_vec,x,useRegul,integrals,sgn,constants,...
                    psiType,useCBIE,useHBIE,dXIdv,dR_xdxi,dR_xdeta,v_1,v_2,alpha);
    for j = 1:n_en
        A_row(sctr_x(j)) = A_row(sctr_x(j)) + R_xScaled(j);
    end     
    A(i,:) = A_row;
    FF(i,:) = getF_eTemp(FF_temp,useNeumanProj,solveForPtot,psiType,useCBIE,useHBIE,useRegul,R_x,sctr_x,x,nx,...
                U,dU,p_inc,dp_inc,dpdn,alpha,integrals,k,constants,sgn);
end
task.varCol{1}.A_K = A;
task.varCol{1}.FF = FF;
% totNoQP
task.varCol{1}.totNoQPnonPolar = totNoQPnonPolar;
task.varCol{1}.totNoQP = totNoQP;
