function F = applyHWBC_ROM(varCol,noVecs)

%coordinateSystem: 0 is cartesian, 1 is cylindrical, 2 is spherical

p_inc = varCol.p_inc;
d_vec = varCol.d_vec;
k = varCol.k;

gluedNodes = varCol.gluedNodes;

noDofs_tot = varCol.noDofs_tot;

weights = varCol.weights;
controlPts = varCol.controlPts;

% rightHandedOrientation = findOrientation(varCol);

elRangeXi = varCol.elRange{1};
elRangeEta = varCol.elRange{2};
patches = varCol.patches;
knotVecs = varCol.knotVecs;
noPatches = varCol.noPatches;

p_xi = varCol.degree(1); % assume p_xi is equal in all patches
p_eta = varCol.degree(2); % assume p_eta is equal in all patches
n_en = (p_xi+1)*(p_eta+1);



% Computes Neumann conditions if the analytic functions is only known at
% the boundary, g_xi0, g_xi1 etc.

F = zeros(noDofs_tot,noVecs);        % external force vector


noParams = 2;
noSurfDofs = 0;
noElemsPatch = zeros(noPatches,1);
noEl = zeros(noPatches,noParams);
for i = 1:noPatches
    n_xi = patches{i}.nurbs.number(1);
    n_eta = patches{i}.nurbs.number(2);
    noSurfDofs = noSurfDofs + n_xi*n_eta;
    for j = 1:noParams
        noEl(i,j) = size(patches{i}.elRange{j},1);
    end
    noElemsPatch(i) = size(patches{i}.elRange{1},1)*size(patches{i}.elRange{2},1);
end

zeta0Nodes = zeros(1,noSurfDofs);
counter = 1;
shiftIdx = 0;
for patch = 1:noPatches
    n_xi = patches{patch}.nurbs.number(1);
    n_eta = patches{patch}.nurbs.number(2);
    n_zeta = patches{patch}.nurbs.number(3);
    for j = 1:n_eta
        for i = 1:n_xi
            zeta0Nodes(counter) = shiftIdx + n_xi*(j-1) + i;
            counter = counter + 1;
        end
    end
    shiftIdx = shiftIdx + n_zeta*n_eta*n_xi;
end
noElems = sum(noElemsPatch);
pIndex = zeros(noElems,1);
element = zeros(noElems,n_en);
index = zeros(noElems,noParams);
e = 1;
maxDof = 0;
jEl = zeros(1,2);
for i = 1:noPatches
    Xi = knotVecs{i}{1};
    Eta = knotVecs{i}{2};
    n_xi = patches{i}.nurbs.number(1);
    n_eta = patches{i}.nurbs.number(2);
    [surfElement, indexXiEta, noElemsXiEta] = generateIGA2DMesh(Xi, Eta, p_xi, p_eta, n_xi, n_eta);
    index(e:e+noElemsPatch(i)-1,:) = indexXiEta + repmat(jEl,noElemsXiEta,1);    
    pIndex(e:e+noElemsXiEta-1) = i;
    element(e:e+noElemsXiEta-1,:) = maxDof + surfElement;
    jEl = jEl + noEl(i,:);
    maxDof = maxDof + n_xi*n_eta;
    e = e + noElemsPatch(i);
end
element2 = element;
% Glue nodes in 2D mesh
for i = 1:length(gluedNodes)
    indices = (zeta0Nodes(element(:)) == gluedNodes{i}(1));
    parentIdx = element(indices);
    element(indices) = parentIdx;
    for j = 2:length(gluedNodes{i})
        indices = (zeta0Nodes(element(:)) == gluedNodes{i}(j));
        element(indices) = parentIdx;
    end
end

Fvalues = zeros(n_en,noElems,noVecs);
indices = zeros(n_en,noElems);

[W2D,Q2D] = gaussianQuadNURBS(p_xi+1,p_eta+1); 

parfor e = 1:noElems
% for e = 1:noElems
    patch = pIndex(e); % New
    Xi = knotVecs{patch}{1}; % New
    Eta = knotVecs{patch}{2}; % New

    idXi = index(e,1);   % the index matrix is made in generateIGA3DMesh
    idEta = index(e,2);

    Xi_e = elRangeXi(idXi,:); % [eta_j,eta_j+1]
    Eta_e = elRangeEta(idEta,:); % [zeta_k,zeta_k+1]

    J_2 = 0.25*(Xi_e(2)-Xi_e(1))*(Eta_e(2)-Eta_e(1));

    sctrXiEta = zeta0Nodes(element(e,:));

    pts = controlPts(sctrXiEta,:);
    wgts = weights(zeta0Nodes(element2(e,:)),:); % New
    
    f_e = zeros(n_en,noVecs);
    for gp = 1:size(W2D,1)
        pt = Q2D(gp,:);
        wt = W2D(gp);

        xi  = parent2ParametricSpace(Xi_e, pt(1));
        eta = parent2ParametricSpace(Eta_e,pt(2));

        [R, dRdxi, dRdeta] = NURBS2DBasis(xi, eta, p_xi, p_eta, Xi, Eta, wgts);
        J = pts'*[dRdxi' dRdeta'];

        
        crossProd = cross(J(:,1), J(:,2)); 


        X = R*pts;
        n = -crossProd/norm(crossProd);
        m = 0:(noVecs-1);
        d_vecX = dot(d_vec,X);
        temp = m/d_vecX;
        temp(1) = 0;
        deriv = -exp(1i*pi*m/2)*dot(d_vec,n)*p_inc(X).*d_vecX.^m.*(temp+1i*k);
        f_e = f_e + R'*deriv*norm(crossProd) * J_2 * wt;  
    end    
    indices(:,e) = sctrXiEta';
    Fvalues(:,e,:) = f_e;
end

for i = 1:noVecs
    F(:,i) = vectorAssembly(Fvalues(:,:,i),indices,noDofs_tot);
end
