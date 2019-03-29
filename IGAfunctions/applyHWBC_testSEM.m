function F = applyHWBC_testSEM(varCol,no_angles)

%coordinateSystem: 0 is cartesian, 1 is cylindrical, 2 is spherical

dp_inc = varCol.dp_inc;

gluedNodes = varCol.gluedNodes;

noDofs_tot = varCol.noDofs_tot;

controlPts = varCol.controlPts;

% rightHandedOrientation = findOrientation(varCol);

patches = varCol.patches;
knotVecs = varCol.knotVecs;
noPatches = varCol.noPatches;

p_xi = varCol.degree(1); % assume p_xi is equal in all patches
p_eta = varCol.degree(2); % assume p_eta is equal in all patches
n_en = (p_xi+1)*(p_eta+1);



% Computes Neumann conditions if the analytic functions is only known at
% the boundary, g_xi0, g_xi1 etc.

F = zeros(noDofs_tot,no_angles);        % external force vector


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

Fvalues = zeros(n_en,noElems,no_angles);
indices = zeros(n_en,noElems);

nurbs = varCol.patches{1}.nurbs;
Nxi = nurbs.number(1);
Neta = nurbs.number(2);
GLL = nurbs.GLL;
tXi = GLL{1};
tEta = GLL{2};
rho = nurbs.rho;
rho_xi = rho{1};
rho_eta = rho{2};
Q2D = [copyVector(tXi.',numel(tEta),1), copyVector(tEta.',numel(tXi),2)];
W2D = copyVector(rho_xi.',numel(tEta),1).*copyVector(rho_eta.',numel(tXi),2);

% parfor e = 1:noElems
for e = 1:noElems
    J_2 = 1;

    sctrXiEta = zeta0Nodes(element(e,:));

    pts = controlPts(sctrXiEta,:);
    
    f_e = zeros(n_en,no_angles);
    for gp = 1:size(W2D,1)
        pt = Q2D(gp,:);
        wt = W2D(gp);

        xi  = pt(1);
        eta = pt(2);
        
%         [R, dRdxi, dRdeta] = NURBS2DBasis(xi, eta, p_xi, p_eta, Xi, Eta, wgts);
        Rxi = zeros(1,Nxi);
        dRdxit = zeros(1,Nxi);
        for i = 1:Nxi
            dRdxit(i) = lagrangePolynomialsDeriv(xi,i,Nxi,tXi);
            Rxi(i) = lagrangePolynomials(xi,i,Nxi,tXi);
        end
        Reta = zeros(1,Neta);
        dRdetat = zeros(1,Neta);
        for j = 1:Neta
            dRdetat(j) = lagrangePolynomialsDeriv(eta,j,Neta,tEta);
            Reta(j) = lagrangePolynomials(eta,j,Neta,tEta);
        end
        R = (copyVector(Rxi,Neta,1).*copyVector(Reta,Nxi,2)).';
        dRdeta = (copyVector(Rxi,Neta,1).*copyVector(dRdetat,Nxi,2)).';
        dRdxi = (copyVector(dRdxit,Neta,1).*copyVector(Reta,Nxi,2)).';
        
        J = pts'*[dRdxi' dRdeta'];

        
        crossProd = cross(J(:,1), J(:,2)); 


        X = R*pts;
        n = -crossProd/norm(crossProd);
                
        deriv = -dp_inc(X,n).';
        f_e = f_e + R'*deriv*norm(crossProd) * J_2 * wt;  
    end    
    indices(:,e) = sctrXiEta';
    Fvalues(:,e,:) = f_e;
end

for alpha_s_Nr = 1:no_angles
    F(:,alpha_s_Nr) = vectorAssembly(Fvalues(:,:,alpha_s_Nr),indices,noDofs_tot);
end
