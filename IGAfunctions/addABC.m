function A = addABC(varCol)

elRangeXi = varCol.elRange{1};
elRangeEta = varCol.elRange{2};
patches = varCol.patches;
knotVecs = varCol.knotVecs;
noPatches = varCol.noPatches;

p_xi = varCol.degree(1); % assume p_xi is equal in all patches
p_eta = varCol.degree(2); % assume p_eta is equal in all patches
n_en = (p_xi+1)*(p_eta+1);

gluedNodes = varCol.gluedNodes;
N = varCol.N;
formulation = varCol.formulation;
A_2 = varCol.A_2;
x_0 = varCol.x_0;
noDofs = varCol.noDofs;

weights = varCol.weights;
controlPts = varCol.controlPts;

k = varCol.k;
Upsilon = varCol.Upsilon;
r_a = varCol.r_a;


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

zeta1Nodes = zeros(1,noSurfDofs);
counter = 1;
shiftIdx = 0;
for patch = 1:noPatches
    n_xi = patches{patch}.nurbs.number(1);
    n_eta = patches{patch}.nurbs.number(2);
    n_zeta = patches{patch}.nurbs.number(3);
    for j = 1:n_eta
        for i = 1:n_xi
            zeta1Nodes(counter) = shiftIdx + (n_eta*n_xi)*(n_zeta-1) + n_xi*(j-1) + i;
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
    indices = (zeta1Nodes(element(:)) == gluedNodes{i}(1));
    parentIdx = element(indices);
    element(indices) = parentIdx;
    for j = 2:length(gluedNodes{i})
        indices = (zeta1Nodes(element(:)) == gluedNodes{i}(j));
        element(indices) = parentIdx;
    end
end
%% Calculate contribution from infinite elements
n_en = (p_xi+1)*(p_eta+1);

spIdxRow = zeros(n_en^2,noElems);
spIdxCol = zeros(n_en^2,noElems);
Avalues = zeros(n_en^2,noElems);

[W2D,Q2D] = gaussianQuadNURBS(p_xi+1,p_eta+1); 
for e = 1:noElems
% parfor e = 1:noElems
    patch = pIndex(e); % New
    Xi = knotVecs{patch}{1}; % New
    Eta = knotVecs{patch}{2}; % New
    idXi = index(e,1); 
    idEta = index(e,2);

    Xi_e = elRangeXi(idXi,:);
    Eta_e = elRangeEta(idEta,:);

    J_2 = 0.25*(Xi_e(2)-Xi_e(1))*(Eta_e(2)-Eta_e(1));

    sctrLocal = element(e,:);
    sctrGlobal = zeta1Nodes(sctrLocal);

    pts = controlPts(sctrGlobal,:);
    wgts = weights(zeta1Nodes(element2(e,:)),:); % New
    
    A_IJ = zeros(n_en);
    
    for gp = 1:size(W2D,1)
        pt = Q2D(gp,:);
        wt = W2D(gp);

        xi  = parent2ParametricSpace(Xi_e, pt(1));
        eta = parent2ParametricSpace(Eta_e,pt(2));
        
        [R, dRdxi, dRdeta] = NURBS2DBasis(xi, eta, p_xi, p_eta, Xi, Eta, wgts);
%         [R, dRdxi, dRdeta, d2Rdxi2, d2Rdeta2] = NURBS2DBasis2(xi, eta, p_xi, p_eta, Xi, Eta, wgts);
        J = pts'*[dRdxi' dRdeta'];
        
        X = R*pts;
        Xt = A_2*(X-x_0)';
        xt = Xt(1);
        yt = Xt(2);
        zt = Xt(3);
        
        [~, theta, ~, c1, c2] = evaluateProlateCoords(xt,yt,zt,Upsilon);
        
        DPDX = dPdX(xt,yt,zt,Upsilon,r_a,c1,c2)*A_2;
        
        J3 = DPDX(2:3,:)*J(:,1:2);
        J_3 = abs(det(J3));
        
        dRdP = J3'\[dRdxi; dRdeta];
        dRdtheta = dRdP(1,:);
        dRdphi = dRdP(2,:);

        switch formulation
            case 'HH'
                switch N
                    case 1
                        A_IJ = A_IJ - (1i*k-1/r_a)*(R'*R)*sqrt(r_a^2-Upsilon^2*cos(theta)^2)*sqrt(r_a^2-Upsilon^2)*sin(theta)*J_3*J_2*wt;   
%                         A_IJ = A_IJ - (1i*k-1/r_a)*(R'*R)*J_3*J_2*wt;   
                    case 2
                        A_IJ = A_IJ - ((1i*k-1/r_a)*(R'*R) ...
                                       - 1/(2*r_a^3)/(1-1i*k*r_a)*(dRdtheta'*dRdtheta + dRdphi'*dRdphi/sin(theta)^2)) ...
                                      *sqrt(r_a^2-Upsilon^2*cos(theta)^2)*sqrt(r_a^2-Upsilon^2)*sin(theta)*J_3*J_2*wt; 
%                         A_IJ = A_IJ + ((1i*k-1/r_a)*(R'*R) ...
%                                        - 1/(2*r_a)/(1-1i*k*r_a)*(R'*dR2dtheta2 + cot(theta)*R'*dRdtheta + R'*dR2dphi2/sin(theta)^2)) ...
%                                       *sqrt(r_a^2-Upsilon^2*cos(theta)^2)*sqrt(r_a^2-Upsilon^2)*sin(theta)*J_3*J_2*wt;  
                end
            case 'BGT'
                switch N
                    case 1
                        A_IJ = A_IJ - (1i*k-1/r_a)*(R'*R)*sqrt(r_a^2-Upsilon^2*cos(theta)^2)*sqrt(r_a^2-Upsilon^2)*sin(theta)*J_3*J_2*wt;  
                    case 2
                        A_IJ = A_IJ - ((1i*k-1/r_a)*(R'*R) ...
                                       - 1/(2*r_a^3)/(1-1i*k*r_a)*(dRdtheta'*dRdtheta + dRdphi'*dRdphi/sin(theta)^2)) ...
                                      *sqrt(r_a^2-Upsilon^2*cos(theta)^2)*sqrt(r_a^2-Upsilon^2)*sin(theta)*J_3*J_2*wt; 
%                         A_IJ = A_IJ + ((1i*k-1/r_a)*(R'*R) ...
%                                        - 1/(2*r_a)/(1-1i*k*r_a)*(R'*dR2dtheta2 + cot(theta)*R'*dRdtheta + R'*dR2dphi2/sin(theta)^2)) ...
%                                       *sqrt(r_a^2-Upsilon^2*cos(theta)^2)*sqrt(r_a^2-Upsilon^2)*sin(theta)*J_3*J_2*wt;  
                end
                
        end
    end   
    
    Avalues(:,e) = reshape(A_IJ,n_en^2,1);
    spIdxRow(:,e) = copyVector(sctrGlobal,n_en,1);
    spIdxCol(:,e) = copyVector(sctrGlobal,n_en,2);
end

spIdxRow = reshape(spIdxRow,numel(spIdxRow),1);
spIdxCol = reshape(spIdxCol,numel(spIdxCol),1);
Avalues = reshape(Avalues,numel(Avalues),1);

[spIdx,~,IuniqueIdx] = unique([spIdxRow, spIdxCol],'rows');
Avalues = accumarray(IuniqueIdx,Avalues);

A = sparse(spIdx(:,1),spIdx(:,2),Avalues,noDofs,noDofs,numel(IuniqueIdx));
