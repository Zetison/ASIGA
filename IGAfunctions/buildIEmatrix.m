function [A, newDofsToRemove, A_1, A_2] = buildIEmatrix(varCol)

elRange = varCol.elRange(1:2);
knotVecs = varCol.knotVecs;

degree = varCol.degree(1:2); % assume degree is equal in all patches

N = varCol.N;
formulation = varCol.formulation;
A_2 = varCol.A_2;
x_0 = varCol.x_0;
noDofs = varCol.noDofs;

weights = varCol.weights;
controlPts = varCol.controlPts;

k = varCol.k(1);
Upsilon = varCol.Upsilon;
r_a = varCol.r_a;

[D,Dt] = generateCoeffMatrix(varCol);
if strcmp(varCol.method,'IENSG')
    noElems = varCol.noElems;
    element = varCol.element;
    element2 = varCol.element2;
    index = varCol.index;
    pIndex = varCol.pIndex;
    n_en = prod(degree+1);
    noSurfDofs = noDofs;
    noDofs = 0;
    nodes = 1:noSurfDofs;
    noDofs_new = noSurfDofs*N;
else
    [nodes, noElems, element, element2, index, pIndex, n_en, noSurfDofs] = meshBoundary(varCol,1);
    noDofs_new = noDofs + noSurfDofs*(N-1);
end


%% Calculate contribution from infinite elements

spIdxRow = zeros(n_en^2,noElems);
spIdxCol = zeros(n_en^2,noElems);
A1values = zeros(n_en^2,noElems);
A2values = zeros(n_en^2,noElems);
A3values = zeros(n_en^2,noElems);
A4values = zeros(n_en^2,noElems);
A5values = zeros(n_en^2,noElems);

[Q, W] = gaussTensorQuad(degree+1);

progressBars = varCol.progressBars;
nProgressStepSize = ceil(noElems/1000);
if progressBars
    ppm = ParforProgMon('Building infinite element matrix: ', noElems, nProgressStepSize);
else
    ppm = NaN;
end

% for e = 1:noElems
parfor e = 1:noElems
	if progressBars && mod(e,nProgressStepSize) == 0
        ppm.increment();
	end
    patch = pIndex(e);
    knots = knotVecs{patch}(1:2);
    Xi_e = zeros(2,2);
    for i = 1:2
        Xi_e(i,:) = elRange{i}(index(e,i),:);
    end
    
    J_2 = prod(Xi_e(:,2)-Xi_e(:,1))/2^2;
    
    xi = parent2ParametricSpace(Xi_e, Q);
    I = findKnotSpans(degree, xi(1,:), knots);
    
    
    sctrLocal = element(e,:);
    sctrGlobal = nodes(sctrLocal);

    pts = controlPts(sctrGlobal,:);
    wgts = weights(nodes(element2(e,:)),:); % New
    
    A1_IJ = zeros(n_en);
    A2_IJ = zeros(n_en);
    A4_IJ = zeros(n_en);
    A5_IJ = zeros(n_en);
    A3_IJ = zeros(n_en);
            
    R = NURBSbasis(I, xi, degree, knots, wgts);
    dXdxi = R{2}*pts;
    dXdeta = R{3}*pts;
    a11 = dXdxi(:,1);
    a21 = dXdxi(:,2);
    a31 = dXdxi(:,3);
    a12 = dXdeta(:,1);
    a22 = dXdeta(:,2);
    a32 = dXdeta(:,3);

    X = R{1}*pts;
    Xt = (X-x_0)*A_2.';

    [~, theta_arr, ~, d1, d2] = evaluateProlateCoords(Xt,Upsilon);
        
    for gp = 1:numel(W)
        theta = theta_arr(gp);
        RR = R{1}(gp,:)'*R{1}(gp,:);
        if sin(theta) < 10*eps
            temp = RR*(a11(gp)*a22(gp)-a12(gp)*a21(gp))/(r_a^2-Upsilon^2);
            temp2 = ( (a11(gp)^2+a21(gp)^2)*(R{3}(gp,:)'*R{3}(gp,:)) ...
                     -(a11(gp)*a12(gp)+a21(gp)*a22(gp))*(R{3}(gp,:)'*R{2}(gp,:) + R{2}(gp,:)'*R{3}(gp,:)) ...
                     +(a12(gp)^2+a22(gp)^2)*(R{2}(gp,:)'*R{2}(gp,:)))/(a11(gp)*a22(gp)-a12(gp)*a21(gp));
        end
        if theta < 10*eps
            A1_IJ = A1_IJ + temp*J_2*W(gp); 
            A2_IJ = A2_IJ + temp2*J_2*W(gp);   
            A3_IJ = A3_IJ + temp*J_2*W(gp); % note that cos(theta)^2 = 1 in this case
        elseif theta > pi-10*eps
            A1_IJ = A1_IJ - temp*J_2*W(gp); 
            A2_IJ = A2_IJ - temp2*J_2*W(gp);   
            A3_IJ = A3_IJ - temp*J_2*W(gp); % note that cos(theta)^2 = 1 in this case
        else
            DPDX = dPdX(Xt(gp,:),Upsilon,r_a,d1(gp),d2(gp))*A_2;
            J = [a11(gp) a12(gp);
                 a21(gp) a22(gp);
                 a31(gp) a32(gp)];
                 
            J3 = DPDX(2:3,:)*J;
            J_3 = abs(det(J3));

            dRdP = J3'\[R{2}(gp,:); R{3}(gp,:)];
            dRdtheta = dRdP(1,:);
            dRdphi = dRdP(2,:);
            
            A1_IJ = A1_IJ + RR*sin(theta)                          	*J_3*J_2*W(gp);  
            A2_IJ = A2_IJ + dRdtheta'*dRdtheta*sin(theta)        	*J_3*J_2*W(gp);  
            A3_IJ = A3_IJ + RR*cos(theta)^2*sin(theta)           	*J_3*J_2*W(gp);  
            A4_IJ = A4_IJ + dRdphi'*dRdphi/sin(theta)              	*J_3*J_2*W(gp);  
            A5_IJ = A5_IJ + dRdphi'*dRdphi*cos(theta)^2/sin(theta) 	*J_3*J_2*W(gp);  
        end
    end   
    
    A1values(:,e) = A1_IJ(:);
    A2values(:,e) = A2_IJ(:);
    A3values(:,e) = A3_IJ(:);
    A4values(:,e) = A4_IJ(:);
    A5values(:,e) = A5_IJ(:);
    spIdxRow(:,e) = copyVector(sctrLocal,n_en,1);
    spIdxCol(:,e) = copyVector(sctrLocal,n_en,2);
end
spIdxRow = reshape(spIdxRow,numel(spIdxRow),1);
spIdxCol = reshape(spIdxCol,numel(spIdxCol),1);
A1values = reshape(A1values,numel(A1values),1);
A2values = reshape(A2values,numel(A2values),1);
A3values = reshape(A3values,numel(A3values),1);
A4values = reshape(A4values,numel(A4values),1);
A5values = reshape(A5values,numel(A5values),1);
[spIdx,~,IuniqueIdx] = unique([spIdxRow, spIdxCol],'rows');
A1values = accumarray(IuniqueIdx,A1values);
A2values = accumarray(IuniqueIdx,A2values);
A3values = accumarray(IuniqueIdx,A3values);
A4values = accumarray(IuniqueIdx,A4values);
A5values = accumarray(IuniqueIdx,A5values);

A1values = sparse(spIdx(:,1),spIdx(:,2),A1values,noSurfDofs,noSurfDofs,numel(IuniqueIdx));
A2values = sparse(spIdx(:,1),spIdx(:,2),A2values,noSurfDofs,noSurfDofs,numel(IuniqueIdx));
A3values = sparse(spIdx(:,1),spIdx(:,2),A3values,noSurfDofs,noSurfDofs,numel(IuniqueIdx));
A4values = sparse(spIdx(:,1),spIdx(:,2),A4values,noSurfDofs,noSurfDofs,numel(IuniqueIdx));
A5values = sparse(spIdx(:,1),spIdx(:,2),A5values,noSurfDofs,noSurfDofs,numel(IuniqueIdx));


%% Evaluate analytic integrals in ``radial'' direction. 
% Note that the last two integrals (B1(end) and B1(end-1),
% B2(end) and B2(end-1)) will be redundant for the cases 'BGC' and 'BGU'
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
ntpmt = nt+mt;
switch formulation
    case 'PGU'
        K1 = -2*varrho2^2*B1(ntpmt) - 1i*varrho2*(ntpmt+2).*B1(ntpmt+1) + ((nt+2)*mt + varrho3^2).*B1(ntpmt+2) ...
             +1i*varrho1*varrho3*(ntpmt+2)*B1(ntpmt+3) - varrho1^2*(nt+2)*mt.*B1(ntpmt+4);
        K2 = B1(ntpmt+2);
        K3 = varrho3^2*B1(ntpmt+2);
        K4 = B2(ntpmt+1);
        K5 = -varrho1^2*B2(ntpmt+3);
    case 'PGC'
        K1_0 = (nt+2)*mt.*B1(ntpmt+2) - varrho1^2*(nt+2)*mt.*B1(ntpmt+4);
        K1_1 = -1i*varrho2*(nt-mt+2).*B1(ntpmt+1) +1i*varrho1*varrho3*(nt-mt+2).*B1(ntpmt+3);
        K1_2 = -varrho3^2.*B1(ntpmt+2);
        K2 = B1(ntpmt+2);
        K3 = varrho3^2*B1(ntpmt+2);
        K4 = B2(ntpmt+1);
        K5 = -varrho1^2*B2(ntpmt+3);
    case 'BGU'
        ntpmt(1) = 3; % To avoid exceeding array bounds (nt,mt=1 is treated separately)
        K1 = -2*varrho2^2*B1(ntpmt-2) - 1i*varrho2*ntpmt.*B1(ntpmt-1) + (nt*mt + varrho3^2).*B1(ntpmt) ...
             +1i*varrho1*varrho3*ntpmt.*B1(ntpmt+1) - varrho1^2*nt*mt.*B1(ntpmt+2);
        K2 = B1(ntpmt);
        K3 = varrho3^2*B1(ntpmt);
        K4 = B2(ntpmt-1);
        K5 = -varrho1^2*B2(ntpmt+1);
        
        K1(1) = -2*1i*varrho2*B1(1) + (1 + varrho3^2)*B1(2) ...
                +2*1i*varrho1*varrho3*B1(3) - varrho1^2*B1(4) - 1i*varrho2*exp(2*1i*varrho2);
        K2(1) = B1(2);
        K3(1) = varrho3^2*B1(2);
        K4(1) = B2(1);
        K5(1) = -varrho1^2*B2(3);
    case 'BGC'
        ntpmt(1) = 3; % To avoid exceeding array bounds (nt,mt=1 is treated separately)
        K1_0 = nt*mt.*B1(ntpmt) - varrho1^2*nt*mt.*B1(ntpmt+2);
        K1_1 = -1i*varrho2*(nt-mt).*B1(ntpmt-1) + 1i*varrho1*varrho3*(nt-mt).*B1(ntpmt+1);
        K1_2 = -varrho3^2.*B1(ntpmt);
        
        
        K2 = B1(ntpmt);
        K3 = varrho3^2*B1(ntpmt);
        K4 = B2(ntpmt-1);
        K5 = -varrho1^2*B2(ntpmt+1);
        
        K1_0(1) = B1(2) - varrho1^2*B1(4);
        K1_1(1) = -1i*varrho2;
        K1_2(1) = -varrho3^2.*B1(2);
        K2(1) = B1(2);
        K3(1) = varrho3^2*B1(2);
        K4(1) = B2(1);
        K5(1) = -varrho1^2*B2(3);
end
switch formulation
    case {'PGC', 'BGC'}
        K1_0 = K1_0*r_a;
        K1_1 = K1_1*r_a;
        K1_2 = K1_2*r_a;
        K2 = K2*r_a;
        K3 = K3*r_a;
        K4 = K4*r_a;
        K5 = K5*r_a;
    case {'PGU', 'BGU'}
        K1 = K1*r_a*exp(-2*1i*varrho2);
        K2 = K2*r_a*exp(-2*1i*varrho2);
        K3 = K3*r_a*exp(-2*1i*varrho2);
        K4 = K4*r_a*exp(-2*1i*varrho2);
        K5 = K5*r_a*exp(-2*1i*varrho2);
end
if varCol.useROM
    A_2 =  kron(Dt*K1_2/k^2*D.',A1values) ...
         + kron(Dt*K3/k^2*D.',A3values);
    A_1 =  kron(Dt*K1_1/k*D.',A1values);
    A =    kron(Dt*K1_0*D.',A1values) ...
         + kron(Dt*K2*D.',A2values) ...
         + kron(Dt*K4*D.',A4values) ...
         + kron(Dt*K5*D.',A5values);
else
    switch formulation
        case {'PGC', 'BGC'}
            K1 = K1_0 + K1_1 + K1_2;
    end
    A  =   kron(Dt*K1*D.',A1values) ...
         + kron(Dt*K2*D.',A2values) ...
         + kron(Dt*K3*D.',A3values) ...
         + kron(Dt*K4*D.',A4values) ...
         + kron(Dt*K5*D.',A5values);
end
  
dofsToRemove = setdiff(1:noSurfDofs,unique(spIdx(:,1)));
noDofsToRemove = numel(dofsToRemove);
newDofsToRemove = zeros(1,noDofsToRemove*N);
i = 1;
for n = 1:N
    newDofsToRemove(i:i+noDofsToRemove-1) = dofsToRemove+(n-1)*noSurfDofs;
    i = i + noDofsToRemove;
end

if noDofs > 0
    [i,j,Avalues] = find(A);    
    [i,j,newDofsToRemove] = modifyIndices(i,j,noSurfDofs,noDofs,nodes,newDofsToRemove);    
    A = sparse(i,j,Avalues,noDofs_new,noDofs_new);
    
    if varCol.useROM
        [i,j,Avalues] = find(A_1);
        [i,j] = modifyIndices(i,j,noSurfDofs,noDofs,nodes);
        A_1 = sparse(i,j,Avalues,noDofs_new,noDofs_new);
        
        [i,j,Avalues] = find(A_2);
        [i,j] = modifyIndices(i,j,noSurfDofs,noDofs,nodes);
        A_2 = sparse(i,j,Avalues,noDofs_new,noDofs_new);
    end
end

function [i,j,newDofsToRemove] = modifyIndices(i,j,noSurfDofs,noDofs,nodes,newDofsToRemove)

indices = i <= noSurfDofs;
indices2 = i > noSurfDofs;
i(indices) = nodes(i(indices));
i(indices2) = i(indices2)+noDofs-noSurfDofs;
indices = j <= noSurfDofs;
indices2 = j > noSurfDofs;
j(indices) = nodes(j(indices));
j(indices2) = j(indices2)+noDofs-noSurfDofs;
if nargout > 2
    indices = newDofsToRemove <= noSurfDofs;
    indices2 = newDofsToRemove > noSurfDofs;
    newDofsToRemove(indices) = nodes(newDofsToRemove(indices));
    newDofsToRemove(indices2) = newDofsToRemove(indices2)+noDofs-noSurfDofs;
end



