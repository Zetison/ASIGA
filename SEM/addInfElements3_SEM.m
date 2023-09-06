function [BB,element,elementInf,zeta1Nodes,zeta0Nodes] = addInfElements3_SEM(task)

patches = task.varCol{1}.patches;
knotVecs = task.varCol{1}.knotVecs;
noPatches = task.varCol{1}.noPatches;

p_xi = task.varCol{1}.degree(1); % assume p_xi is equal in all patches
p_eta = task.varCol{1}.degree(2); % assume p_eta is equal in all patches
n_en = (p_xi+1)*(p_eta+1);

gluedNodes = task.varCol{1}.gluedNodes;
N = task.iem.N;
formulation = task.misc.formulation;

k = task.misc.omega/task.varCol{1}.c_f;
Upsilon = task.iem.Upsilon;
r_a = task.misc.r_a;


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
elementInf = zeros(noElems,n_en);
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
    elementInf(e:e+noElemsXiEta-1,:) = maxDof + surfElement;
    jEl = jEl + noEl(i,:);
    maxDof = maxDof + n_xi*n_eta;
    e = e + noElemsPatch(i);
end
% Glue nodes in 2D mesh
for i = 1:length(gluedNodes)
    indices = (zeta1Nodes(elementInf(:)) == gluedNodes{i}(1));
    parentIdx = elementInf(indices);
    elementInf(indices) = parentIdx;
    for j = 2:length(gluedNodes{i})
        indices = (zeta1Nodes(elementInf(:)) == gluedNodes{i}(j));
        elementInf(indices) = parentIdx;
    end
end


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


%% Evaluate analytic integrals in ``radial'' direction. 
% Note that the last two integrals (I1(end) and I1(end-1),
% I2(end) and I2(end-1)) will be redundant for the cases 'BGC' and 'BGU'
B1 = zeros(2*N+4,1);
B2 = zeros(2*N+3,1);
varrho1 = Upsilon/r_a;
varrho2 = k*r_a;
varrho3 = k*Upsilon;
for n = 1:2*N+4
    B1(n) = radialIntegral3(n, varrho1, varrho2, formulation, 1);
    if n < 2*N+4
        B2(n) = radialIntegral3(n, varrho1, varrho2, formulation, 2);
    end
end
nt = (1:N).';
mt = 1:N;
BB = zeros(N,N,5);
switch formulation
    case 'PGU'
        BB(:,:,1) = -2*varrho2^2*B1(nt+mt) - 1i*varrho2*(nt+mt+2).*B1(nt+mt+1) + ((nt+2)*mt + varrho3^2).*B1(nt+mt+2) ...
                      +1i*varrho1*varrho3*(nt+mt+2).*B1(nt+mt+3) - varrho1^2*((nt+2)*mt).*B1(nt+mt+4);
        BB(:,:,2) = B1(nt+mt+2);
        BB(:,:,3) = varrho3^2*B1(nt+mt+2);
        BB(:,:,4) = B2(nt+mt+1);
        BB(:,:,5) = -varrho1^2*B2(nt+mt+3);
        BB = BB*exp(-2*1i*varrho2)*r_a;
    case 'PGC'
        BB(:,:,1) = - 1i*varrho2*(nt-mt+2).*B1(nt+mt+1) + ((nt+2)*mt - varrho3^2).*B1(nt+mt+2) ...
                      +1i*varrho1*varrho3*(nt-mt+2).*B1(nt+mt+3) - varrho1^2*((nt+2)*mt).*B1(nt+mt+4);
        BB(:,:,2) = B1(nt+mt+2);
        BB(:,:,3) = varrho3^2*B1(nt+mt+2);
        BB(:,:,4) = B2(nt+mt+1);
        BB(:,:,5) = - varrho1^2*B2(nt+mt+3);
        BB = BB*r_a;
    otherwise
        for nt = 1:N
            for mt = 1:N
                switch formulation
                    case 'BGU'
                        if mt+nt == 2
                            BB(nt,mt,1)  = -2*1i*varrho2*B1(1) + (1 + varrho3^2)*B1(2) ...
                                          +2*1i*varrho1*varrho3*B1(3) - varrho1^2*B1(4) - 1i*varrho2*exp(2*1i*varrho2);
                            BB(nt,mt,2)  = B1(2);
                            BB(nt,mt,3)  = varrho3^2*B1(2);
                            BB(nt,mt,4)  = B2(1);
                            BB(nt,mt,5)  = -varrho1^2*B2(3);
                        else
                            BB(nt,mt,1)  = -2*varrho2^2*B1(nt+mt-2) - 1i*varrho2*(nt+mt)*B1(nt+mt-1) + (nt*mt + varrho3^2)*B1(nt+mt) ...
                                          +1i*varrho1*varrho3*(nt+mt)*B1(nt+mt+1) - varrho1^2*nt*mt*B1(nt+mt+2);
                            BB(nt,mt,2)  = B1(nt+mt);
                            BB(nt,mt,3)  = varrho3^2*B1(nt+mt);
                            BB(nt,mt,4)  = B2(nt+mt-1);
                            BB(nt,mt,5)  = -varrho1^2*B2(nt+mt+1);
                        end
                        BB(nt,mt,:) = BB(nt,mt,:)*exp(-2*1i*varrho2)*r_a;
                    case 'BGC'
                        if mt+nt == 2
                            BB(nt,mt,1)  = (1 - varrho3^2)*B1(2) - varrho1^2*B1(4) - 1i*varrho2;
                            BB(nt,mt,2)  = B1(2);
                            BB(nt,mt,3)  = varrho3^2*B1(2);
                            BB(nt,mt,4)  = B2(1);
                            BB(nt,mt,5)  = -varrho1^2*B2(3);
                        else
                            BB(nt,mt,1)  = -1i*varrho2*(nt-mt)*B1(nt+mt-1) + (nt*mt - varrho3^2)*B1(nt+mt) ...
                                          +1i*varrho1*varrho3*(nt-mt)*B1(nt+mt+1) - varrho1^2*nt*mt*B1(nt+mt+2);
                            BB(nt,mt,2)  = B1(nt+mt);
                            BB(nt,mt,3)  = varrho3^2*B1(nt+mt);
                            BB(nt,mt,4)  = B2(nt+mt-1);
                            BB(nt,mt,5)  = -varrho1^2*B2(nt+mt+1);
                        end
                        BB(nt,mt,:) = BB(nt,mt,:)*r_a;
                end
            end
        end
end


     

function [element, index, noElems] = generateIGA2DMesh(Xi, Eta, p, q, n, m)
% warning('This function is depricated. Use generateIGAmesh instead')

uniqueXiKnots   = unique(Xi);
uniqueEtaKnots   = unique(Eta);

noElemsXi       = length(uniqueXiKnots)-1; % # of elements xi dir.
noElemsEta       = length(uniqueEtaKnots)-1; % # of elements eta dir.

chan  = zeros(m,n);

count = 1;
for j = 1:m
    for i = 1:n
        chan(i,j) = count;
        count     = count + 1;
    end
end


% determine our element ranges and the corresponding
% knot indices along each direction

[elRangeXi,elConnXi] = buildConnectivity(p,Xi,noElemsXi);
[elRangeEta,elConnEta] = buildConnectivity(q,Eta,noElemsEta);

noElems = noElemsXi * noElemsEta;
element = zeros(noElems,(p+1)*(q+1));

e = 1;
for jj = 1:noElemsEta
    vConn = elConnEta(jj,:);
    for ii = 1:noElemsXi
        c = 1;
        uConn = elConnXi(ii,:);

        for j = 1:length(vConn)
            for i = 1:length(uConn)
                element(e,c) = chan(uConn(i), vConn(j));
                c = c + 1;
            end
        end
        e = e + 1;
    end
end

index = zeros(noElems,2);

count = 1;
for j = 1:size(elRangeEta,1)
    for i = 1:size(elRangeXi,1)
        index(count,1) = i;
        index(count,2) = j;
        count = count + 1;
    end
end









