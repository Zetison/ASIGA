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
else
    no_angles = 1;
end
solveForPtot = varCol.solveForPtot;
if solveForPtot
    p_inc = varCol.p_inc;
    dp_inc = varCol.dp_inc;
    dpdn = @(x,n) 0;
else
    p_inc = NaN;
    dp_inc = NaN;
    if SHBC
        dpdn = @(x,n) -varCol.dp_inc(x,n);
    else
        dpdn = varCol.dpdn;
    end
end

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
    switch model
        case 'S1'
            Xarr = linspace(-0.5,1,3)/4;
            Yarr = linspace(-0.5,1,3)/4;
            Zarr = linspace(-0.5,1,3)/4;
            [X,Y,Z] = ndgrid(Xarr,Yarr,Zarr);
            cp = [X(:),Y(:), Z(:)];
        case {'BCA','BCA_P'}
            alpha = NaN;
            setBCParameters
            h = g2*tan(alpha/2); % = g2/tan((pi-alpha)/2)
            x2 = g3*tan(alpha);
            h_tail = -b+h+x2;
            xi = [a/4               -L-g3;
                  x_s-l_ls*0.1      x_s-l_ls*0.4;
                  x_d-l_ld*0.1      x_d-l_ld*0.4;
                  x_d-l_ld*0.1      x_d-l_ld*0.4;
                  x_m-l_lm*0.1      x_m-l_lm*0.4;
                  x_m-l_lm*0.1      x_m-l_lm*0.4;
                  x_m-l_lm*0.1      x_m-l_lm*0.4;
                  x_m-l_lm*0.1      x_m-l_lm*0.4;
                  x_m-l_lm*0.9      x_m-l_lm-0.2*(g3+g2+L+x_m-l_lm)];
            yi = [-b/4              b/4;
                  -s/4              s/4;
                  s                 s+0.5*h_d;
                  -s                -s-0.5*h_d;
                  -0.1*b_lm         0.1*b_lm;
                  -0.3*h_m          -0.8*h_m;
                  -0.1*b_lm         0.1*b_lm;
                  0.3*h_m           0.8*h_m;
                  -h_tail*0.5       h_tail*0.5];
            zi = [-b/4              c/4;
                  c              	c+0.5*h_s;
                  c-b_ld/2-0.4*b_ld/2 c-b_ld/2+0.4*b_ld/2;
                  c-b_ld/2-0.4*b_ld/2 c-b_ld/2+0.4*b_ld/2;
                  0.3*h_m           0.8*h_m;
                  -0.1*b_lm         0.1*b_lm;
                  -0.3*h_m          -0.8*h_m;
                  -0.1*b_lm         0.1*b_lm;
                  -h_tail*0.5       h_tail*0.5];
                  
            cp = [];
            npts = 3;
            for i = 1:size(xi,1)
                Xarr = linspace(xi(i,1),xi(i,2),npts);
                Yarr = linspace(yi(i,1),yi(i,2),npts);
                Zarr = linspace(zi(i,1),zi(i,2),npts);
                [X,Y,Z] = ndgrid(Xarr,Yarr,Zarr);
                cp = [cp; X(:),Y(:), Z(:)];
            end
        otherwise
            error('Not implemented')
    end            
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_GP = 0;
pD = plotBEMGeometry(patches,plot_GP,10,0,0.8);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for i = 1:n_cp
%     x = cp(i,:);
%     if plot_GP
%         pD = plotGP(pD,x,'red');
%     end
% end
totNoQP = 0;
% for i = 1:n_cp
parfor i = 1:n_cp
    x = cp(i,:);
    A_row = complex(zeros(1, noDofs));
    [constants, integrals] = initializeBIE(psiType,useRegul,x,NaN,k,model);
    
    FF_temp = complex(zeros(1, no_angles));
    
    for e_y = 1:noElems  
        [BIE, integrals, FF_temp, sctr_y, noGp] = getBEMquadPts(e_y,Q2D_2,W2D_2,Q,W,integrals,FF_temp,...
                useEnrichedBfuns,k,d_vec,useNeumanProj,solveForPtot,useCBIE,useHBIE,dpdn,U,...
                x,NaN,NaN,NaN,NaN,constants,psiType,useRegul,...
                p_xi, p_eta,pIndex,knotVecs,index,elRangeXi,elRangeEta,element,element2,controlPts,weights,...
                patches,Eps,diagsMax,centerPts,agpBEM,quadMethodBEM);
        for j = 1:n_en
            A_row(sctr_y(j)) = A_row(sctr_y(j)) + BIE(j);
        end
        totNoQP = totNoQP + noGp;
    end       
    A(i,:) = A_row;
    FF(i,:) = getF_eTemp(FF_temp,useNeumanProj,solveForPtot,psiType,useCBIE,useHBIE,useRegul,NaN,NaN,x,NaN,...
                U,dU,p_inc,dp_inc,dpdn,NaN,NaN,NaN,NaN,sgn);
end

% totNoQP
varCol.totNoQP = varCol.totNoQP + totNoQP;
