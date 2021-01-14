function plotGalerkinResidual(varCol,Uc)

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
colMethod = varCol.colMethod;

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


p_inc = varCol.p_inc;
Eps = 1e4*eps;

patches = varCol.patches;
[~, ~, diagsMax] = findMaxElementDiameter(patches);
centerPts = findCenterPoints(patches);

noPatches = varCol.noPatches;
x_grev = [];
x_cg = [];
x_gl = [];
for patch = 1:noPatches
    nurbs = patches{patch}.nurbs;
    Xi = nurbs.knots{1};
    Eta = nurbs.knots{2};

    xi1D = splinesGL(Xi,p_xi);
    eta1D = splinesGL(Eta,p_eta);    
    xi = copyVector(xi1D,numel(eta1D),1);
    eta = copyVector(eta1D,numel(xi1D),2);    
    x = evaluateNURBS_2ndDeriv(nurbs, [xi,eta]);
    x_gl = [x_gl; x];
    
    xi1D = aveknt(Xi,p_xi+1).';
    eta1D = aveknt(Eta,p_eta+1).';    
    xi = copyVector(xi1D,numel(eta1D),1);
    eta = copyVector(eta1D,numel(xi1D),2);      
    x = evaluateNURBS_2ndDeriv(nurbs, [xi,eta]);
    x_grev = [x_grev; x];

    xi1D = CauchyGalerkin(p_xi, numel(xi1D), Xi);
    eta1D = CauchyGalerkin(p_eta, numel(eta1D), Eta);
    xi = copyVector(xi1D,numel(eta1D),1);
    eta = copyVector(eta1D,numel(xi1D),2);      
    x = evaluateNURBS_2ndDeriv(nurbs, [xi,eta]);
    x_cg = [x_cg; x];
end

pD.plotGP = true;
pD.plotPointsAsSpheres = true;
pD.pointsRadius = 12e-3;
pD.noSpherePts = 40;
pD.lineColor = 'blue';
pD.lineStyle = '-';
pD.patches = patches;
hold on
pD.h = gca;
pD = plotGP(pD,x_grev,'red');
pD = plotGP(pD,x_gl,'blue');
if mod(p_xi,2)
    pD = plotGP(pD,x_cg,'green');
end
drawnow
keyboard


eNeighbour = NaN; % to avoid transparency "bug"
createElementTopology
% legend show
% set(0,'DefaultLegendAutoUpdate','off')
[Q2D_2,W2D_2,Q,W] = getBEMquadData(p_xi,p_eta,extraGP,extraGPBEM,quadMethodBEM);

for e_x = 1:noElems  
    patch = pIndex(e_x); % New
    Xi = knotVecs{patch}{1}; % New
    Eta = knotVecs{patch}{2}; % New

    idXi = index(e_x,1);
    idEta = index(e_x,2);

    Xi_e_x = elRangeXi(idXi,:);
    Eta_e_x = elRangeEta(idEta,:);
    noElemsXi = numel(unique(Xi))-1;

    sctr = element(e_x,:);
    pts = controlPts(sctr,:);
    wgts = weights(element2(e_x,:)); % New      
    U_sctr = Uc{1}(sctr,:); % New      
    npts = 512/noElemsXi;
%     npts = 256/noElemsXi;
%     npts = 32/noElemsXi;

    xi_x_arr = copyVector(linspace(Xi_e_x(1)+Eps,Xi_e_x(2)-Eps,npts),npts,1);
    eta_x_arr = copyVector(linspace(Eta_e_x(1)+Eps,Eta_e_x(2)-Eps,npts),npts,2);
    [x,u_x] = numericalSolEvalVectorized(xi_x_arr,eta_x_arr,p_xi,p_eta,Xi,Eta,wgts,pts,U_sctr);
    xi_x_arr = copyVector(linspace(Xi_e_x(1),Xi_e_x(2),npts),npts,1);
    eta_x_arr = copyVector(linspace(Eta_e_x(1),Eta_e_x(2),npts),npts,2);
    if 1
%         u_x = u_x + p_inc(x);
        residual = p_inc(x) - u_x;
    else
        u_x = u_x - p_inc(x);
        residual = -u_x;
    end
    parfor i = 1:numel(eta_x_arr)
        xi_x = xi_x_arr(i);
        eta_x = eta_x_arr(i);

%         [constants, integrals] = initializeBIE(psiType,useRegul,x,NaN,k,model);
        [adjacentElements, xi_x_tArr,eta_x_tArr] = getAdjacentElements(e_x,xi_x,eta_x,Xi_e_x,Eta_e_x,eNeighbour,Eps);
        residual_temp = complex(0);
        integrals_temp = complex(0);
%         for e_y = 1:noElems  
        for e_y = 1:noElems  
            [BIE, integrals_temp2, ~, sctr_y] = getBEMquadPts(e_y,Q2D_2,W2D_2,Q,W,complex(0),NaN,...
                    useEnrichedBfuns,k,d_vec,0,SHBC,useCBIE,useHBIE,dpdn,NaN,...
                    x(i,:),NaN(1,3),xi_x_tArr,eta_x_tArr,adjacentElements,NaN,psiType,useRegul,...
                    p_xi, p_eta,pIndex,knotVecs,index,elRangeXi,elRangeEta,element,element2,controlPts,weights,...
                    patches,Eps,diagsMax,centerPts,agpBEM,quadMethodBEM);
            integrals_temp = integrals_temp + integrals_temp2;
            U_sctr = Uc{1}(sctr_y,:); % New      
            residual_temp = residual_temp + BIE*U_sctr;
        end
        residual(i) = residual(i) + residual_temp - integrals_temp*u_x(i);
    end
    X = reshape(x(:,1),npts,npts);
    Y = reshape(x(:,2),npts,npts);
    Z = reshape(x(:,3),npts,npts);
    C = reshape(residual,npts,npts);
%     C = reshape(u_x-(analytic(x)+p_inc(x)),npts,npts);
%     C = reshape(u_x,npts,npts);
%     C = reshape(analytic(x),npts,npts);
%     surf(X,Y,Z, abs(C), 'EdgeColor','none','LineStyle','none')
    max_res = max(max(abs(C)));
    hold on
    surf(X,Y,Z, log10(abs(C)), 'EdgeColor','none','LineStyle','none')
    colorbar 
    axis off
%     caxis([-10,0])
    caxis([-7,-3])
    axis equal
    
    set(gca, 'Color', 'none');
%     title(['Fluid 3D NURBS geometry. Mesh ' num2str(M)])
%     view(-70,30)
%     view(120,10)
    view(18,10)
%     view(0,0)
    camproj('perspective')
    ax = gca;               % get the current axis
    ax.Clipping = 'off';    % turn clipping off
    lineColor = 'black';
    lineWidth = 0.5;
    plot3(X(:,1),Y(:,1),Z(:,1),'color',lineColor,'LineWidth',lineWidth)
    plot3(X(:,end),Y(:,end),Z(:,end),'color',lineColor,'LineWidth',lineWidth)
    plot3(X(1,:),Y(1,:),Z(1,:),'color',lineColor,'LineWidth',lineWidth)
    plot3(X(end,:),Y(end,:),Z(end,:),'color',lineColor,'LineWidth',lineWidth)
    drawnow
end