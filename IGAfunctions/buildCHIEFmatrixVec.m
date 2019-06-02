function [A, FF, varCol] = buildCHIEFmatrixVec(varCol)

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
noDofs = varCol.noDofs;

extraGP = varCol.extraGP;
extraGPBEM = varCol.extraGPBEM;
agpBEM = varCol.agpBEM;
exteriorProblem = varCol.exteriorProblem;
model = varCol.model;

quadMethodBEM = varCol.quadMethodBEM;

Eps = 10*eps;

k = varCol.k;

useCBIE = true;
useHBIE = false;
useRegul = false;
psiType = NaN;

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
patches = varCol.patches;

[~, ~, diagsMax] = findMaxElementDiameter(patches);
centerPts = findCenterPoints(patches);

%% Create source points:
if exteriorProblem
    Xarr = linspace(-0.5,1,3)/4;
    Yarr = linspace(-0.5,1,3)/4;
    Zarr = linspace(-0.5,1,3)/4;
    [X,Y,Z] = ndgrid(Xarr,Yarr,Zarr);
    cp = [X(:),Y(:), Z(:)];
%     cp = [0.49,0,0;
%           0,0,0];
else
    error('Not implemented')
end
n_cp = size(cp,1);

useNeumanProj = varCol.useNeumanProj;
if useNeumanProj
    [U,dU] = projectBC(varCol,SHBC,useCBIE,useHBIE);
else
    U = NaN;
    dU = NaN;
end
n_en = (p_xi+1)*(p_eta+1);

[Q2D_2,W2D_2,Q,W] = getBEMquadData(p_xi,p_eta,extraGP,extraGPBEM,quadMethodBEM);

A = complex(zeros(n_cp, noDofs));
FF = complex(zeros(n_cp, no_angles));

totNoQP = 0;
% for i = 1:n_cp
parfor i = 1:n_cp
    x = cp(i,:);
    A_row = complex(zeros(1, noDofs));
    [constants, integrals] = initializeBIE(psiType,useRegul,x,NaN,k,model);
    
    FF_temp = complex(zeros(1, no_angles));
    
    for e_y = 1:noElems  
        [BIE, integrals, FF_temp, sctr_y, noGp] = getBEMquadPts(e_y,Q2D_2,W2D_2,Q,W,integrals,FF_temp,...
                useEnrichedBfuns,k,d_vec,useNeumanProj,SHBC,useCBIE,useHBIE,dpdn,U,...
                x,NaN,NaN,NaN,NaN,NaN,NaN,constants,psiType,useRegul,...
                p_xi, p_eta,pIndex,knotVecs,index,elRangeXi,elRangeEta,element,element2,controlPts,weights,...
                patches,Eps,diagsMax,centerPts,agpBEM,quadMethodBEM);
        for j = 1:n_en
            A_row(sctr_y(j)) = A_row(sctr_y(j)) + BIE(j);
        end
        totNoQP = totNoQP + noGp;
    end       
    A(i,:) = A_row;
    FF(i,:) = getF_eTemp(FF_temp,useNeumanProj,SHBC,psiType,useCBIE,useHBIE,useRegul,NaN,NaN,x,NaN,...
                U,dU,p_inc,dp_inc,dpdn,NaN,NaN,NaN,NaN,sgn);
end

% totNoQP
varCol.totNoQP = varCol.totNoQP + totNoQP;
