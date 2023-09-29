function task = buildIEmatrix(task)

elRange = task.varCol{1}.elRange(1:2);
knotVecs = task.varCol{1}.knotVecs;

degree = task.varCol{1}.degree(1:2); % assume degree is equal in all patches
Ntot = task.iem.N;
r_a = task.misc.r_a;
if isnan(r_a)
    error('misc.r_a must be set, either in the mesh script or in the task script')
end
IEbasis = task.iem.IEbasis;
p_ie = task.iem.p_ie;
IElocSup = task.iem.IElocSup;
[D,Dt,y_n] = generateCoeffMatrix(task);
if IElocSup
    N = p_ie;
    if Ntot-2*p_ie < 0
        warning('Too low input N for local IE basis functions. Using N = 2*p_ie instead')
        Ntot = 2*p_ie;
    end
    noElems = round(Ntot/p_ie);
    if strcmp(IEbasis,'Lagrange')
        ie_Zeta = [0, reshape(repmat((0:noElems-1)/(noElems-1),p_ie,1),1,noElems*p_ie), 1];
    else
        ie_Zeta = [zeros(1,p_ie+1), linspace2(0,1,Ntot-2*p_ie), ones(1,p_ie+1)];
    end
    rho = r_a*y_n(1:end-p_ie+1);
else
    N = Ntot;
    rho = NaN;
    p_ie = N;
end
formulation = task.misc.formulation;
A_2 = task.iem.A_2;
x_0 = task.iem.x_0;
noDofs = task.varCol{1}.noDofs;

weights = task.varCol{1}.weights;
controlPts = task.varCol{1}.controlPts;

omega = task.misc.omega;
if strcmp(formulation(end),'U') % unconjugated formulation
    k = omega/task.varCol{1}.c_f;
else
    k = 1/task.varCol{1}.c_f; % multiply with task.misc.omega later (in order to facilitate storage of frequency independent matrices)
end
Upsilon = task.iem.Upsilon;
if task.varCol{1}.boundaryMethod
    noElems = task.varCol{1}.noElems;
    element = task.varCol{1}.element;
    element2 = task.varCol{1}.element2;
    index = task.varCol{1}.index;
    pIndex = task.varCol{1}.pIndex;
    n_en = prod(degree+1);
    noSurfDofs = noDofs;
    noDofs = 0;
    nodes = 1:noSurfDofs;
    noDofs_new = noSurfDofs*N;
else
    varColBdry = meshBoundary(task.varCol{1},'Gamma_a');
    
    nodes = varColBdry.nodes;
    noElems = varColBdry.noElems;
    element = varColBdry.element;
    knotVecs = varColBdry.knotVecs;
    elRange = varColBdry.elRange;
    element2 = varColBdry.element2;
    index = varColBdry.index;
    pIndex = varColBdry.pIndex;
    n_en = varColBdry.n_en;
    noSurfDofs = varColBdry.noSurfDofs;
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

extraGP = task.misc.extraGP;
[Q, W] = gaussTensorQuad(degree+1+extraGP(1:2));

progressBars = task.misc.progressBars;
nProgressStepSize = ceil(noElems/1000);
if progressBars
    try
        ppm = ParforProgMon('Building infinite element matrix: ', noElems, nProgressStepSize);
    catch
        progressBars = false;
        ppm = NaN;
    end
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

if IElocSup
    s = task.iem.s_ie;
    x = @(zeta) 1 + zeta.^s*(1/Ntot - 1);
    zeta_a = ((Ntot-p_ie)/(Ntot-1))^(1/s);
    nurbs = parmFunc(ie_Zeta,p_ie,@(xi) x(xi*zeta_a));
    
    varCol_ie.nurbs = nurbs;
    varCol_ie.dimension = 1; 
    varCol_ie = findDofsToRemove(generateIGAmesh(convertNURBS(varCol_ie)),1);
    r_b = r_a./varCol_ie.nurbs{1}.coeffs(1,end);
    if strcmp(IEbasis,'Lagrange')
        r_b = rho(end);
    end
    ra = r_b;
else
    ra = r_a;
end

if 1
    B1 = zeros(2*N+4,1);
    B2 = zeros(2*N+3,1);
    varrho1 = Upsilon/ra;
    varrho2 = k*ra;
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
                 +1i*varrho1*varrho3*(ntpmt+2).*B1(ntpmt+3) - varrho1^2*(nt+2)*mt.*B1(ntpmt+4);
            K2 = B1(ntpmt+2);
            K3 = varrho3^2*B1(ntpmt+2);
            K4 = B2(ntpmt+1);
            K5 = -varrho1^2*B2(ntpmt+3);
        case 'PGC'
            K1 = (nt+2)*mt.*B1(ntpmt+2) - varrho1^2*(nt+2)*mt.*B1(ntpmt+4);
            K1_1 = -1i*varrho2*(nt-mt+2).*B1(ntpmt+1) + 1i*varrho1*varrho3*(nt-mt+2).*B1(ntpmt+3);
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
            K1 = nt*mt.*B1(ntpmt) - varrho1^2*nt*mt.*B1(ntpmt+2);
            K1_1 = -1i*varrho2*(nt-mt).*B1(ntpmt-1) + 1i*varrho1*varrho3*(nt-mt).*B1(ntpmt+1);
            K1_2 = -varrho3^2.*B1(ntpmt);


            K2 = B1(ntpmt);
            K3 = varrho3^2*B1(ntpmt);
            K4 = B2(ntpmt-1);
            K5 = -varrho1^2*B2(ntpmt+1);

            K1(1) = B1(2) - varrho1^2*B1(4);
            K1_1(1) = -1i*varrho2;
            K1_2(1) = -varrho3^2.*B1(2);
            K2(1) = B1(2);
            K3(1) = varrho3^2*B1(2);
            K4(1) = B2(1);
            K5(1) = -varrho1^2*B2(3);
        case 'WBGU'
            K1 = -2*varrho2^2*B1(ntpmt) - 1i*varrho2*ntpmt.*B1(ntpmt+1) + (nt*mt + varrho3^2).*B1(ntpmt+2) ...
                 +1i*varrho1*varrho3*ntpmt.*B1(ntpmt+3) - varrho1^2*nt*mt.*B1(ntpmt+4);
            K2 = B1(ntpmt+2);
            K3 = varrho3^2*B1(ntpmt+2);
            K4 = B2(ntpmt+1);
            K5 = -varrho1^2*B2(ntpmt+3);
        case 'WBGC'
            K1 = nt*mt.*B1(ntpmt+2) - varrho1^2*nt*mt.*B1(ntpmt+4);
            K1_1 = -1i*varrho2*(nt-mt).*B1(ntpmt+1) + 1i*varrho1*varrho3*(nt-mt).*B1(ntpmt+3);
            K1_2 = -varrho3^2.*B1(ntpmt+2);
            K2 = B1(ntpmt+2);
            K3 = varrho3^2*B1(ntpmt+2);
            K4 = B2(ntpmt+1);
            K5 = -varrho1^2*B2(ntpmt+3);
    end
    K1 = Dt*K1*D.'*ra;
    K2 = Dt*K2*D.'*ra;
    K3 = Dt*K3*D.'*ra;
    K4 = Dt*K4*D.'*ra;
    K5 = Dt*K5*D.'*ra;
    switch formulation
        case {'PGC', 'BGC', 'WBGC'}
            K1_1 = Dt*K1_1*D.'*ra;
            K1_2 = Dt*K1_2*D.'*ra;
    end
else
    rho2 = r_a*y_n(end-p_ie+1:end);
%     [Q, W] = gaussTensorQuad(50,'Laguerre');
    [Q, W] = gaussTensorQuad(10);
    varrho3 = k*Upsilon;
    if strcmp(IEbasis,'Lagrange')
        rho_1 = rho2(1);
        r = 2*rho_1./(1-Q);
        J_1 = 2*rho_1./(1-Q).^2;
        x = 1./r;
        x_i = 1./rho2;
        [L,dLdx] = lagrangePolynomials(x,1:p_ie,p_ie,x_i);
        L = L.*x/x_i(1);
        dLdx = dLdx.*x/x_i(1) + L/x_i(1);
        R = cell(1,2);
        R{1} = L;
        dxdr = -1./r.^2;
        dRdr = dLdx.*dxdr;
    else
        zeta = parent2ParametricSpace(Xi_e, Q);
        I = findKnotSpans(degree, zeta(1,:), knots);
        R = NURBSbasis(I, zeta, degree, knots, wgts);

        x = R{1}*pts;
        r = r_a./x;
        dxdzeta = R{2}*pts;
        J_1 = -r_a./x.^2.*dxdzeta; % = drdzeta
        dRdr = R{2}./J_1;
    end
    n_en = p_ie;
    K1 = zeros(n_en);
    K1_1 = zeros(n_en);
    K1_2 = zeros(n_en);
    K2 = zeros(n_en);
    K3 = zeros(n_en);
    K4 = zeros(n_en);
    K5 = zeros(n_en);
    for i = 1:numel(W)
        RR = R{1}(i,:).'*R{1}(i,:);
        RdRdr = R{1}(i,:).'*dRdr(i,:);
        fact = abs(J_1(i)) * W(i);
        rUps = r(i)^2-Upsilon^2;
        kr = k*r(i);
        switch formulation
            case {'WBGC', 'WBGU'}
                fact = fact * (r_a/r(i))^2;
        end
        switch formulation
            case {'PGC', 'BGC'}
                K1   = K1   + rUps*dRdr(i,:).'*dRdr(i,:) * fact;  
                K1_1 = K1_1 + rUps*1i*k*(RdRdr - RdRdr.') * fact;
                K1_2 = K1_2 - varrho3^2*RR                * fact;   
            case {'PGU', 'BGU'}
                fact = fact * exp(2*1i*kr);
                K1 = K1 + (  rUps*dRdr(i,:).'*dRdr(i,:) ... 
                           + rUps*1i*k*(RdRdr + RdRdr.') ...
                           - (2*kr^2-varrho3^2)*RR     ) * fact;  
        end
        K2 = K2 +                  RR * fact;  
        K3 = K3 + varrho3^2      * RR * fact;  
        K4 = K4 + r(i)^2/rUps    * RR * fact;   
        K5 = K5 - Upsilon^2/rUps * RR * fact;  
    end    
end

if IElocSup      
    %% Extract all needed data from options and varCol_ie
    degree = varCol_ie.degree;
    knotVecs = varCol_ie.knotVecs;
    index = varCol_ie.index;
    noElems = varCol_ie.noElems;
    elRange = varCol_ie.elRange;
    element = varCol_ie.element;
    element2 = varCol_ie.element2;
    weights = varCol_ie.weights;
    controlPts = varCol_ie.controlPts;
    pIndex = varCol_ie.pIndex;
    noDofs_ie = varCol_ie.noDofs;
    
    %% Preallocation and initiallizations
    n_en = prod(degree+1);
    sizeKe = n_en^2;
    spIdxRow = zeros(sizeKe,noElems);
    spIdxCol = zeros(sizeKe,noElems);
    K1values = zeros(sizeKe,noElems); 
    K1_1values = zeros(sizeKe,noElems); 
    K1_2values = zeros(sizeKe,noElems); 
    K2values = zeros(sizeKe,noElems); 
    K3values = zeros(sizeKe,noElems); 
    K4values = zeros(sizeKe,noElems); 
    K5values = zeros(sizeKe,noElems); 

%     [Q, W] = gaussTensorQuad(degree+1+extraGP(1));
    [Q, W] = gaussTensorQuad(50);
    %% Build global matrices
    
%     for e = 1:noElems
    parfor e = 1:noElems
        patch = pIndex(e);
        knots = knotVecs{patch};
        Xi_e = elRange{1}(index(e,1),:);
        
        sctr = element(e,:);

        if strcmp(IEbasis,'Lagrange')
            rho_1 = rho(1+(e-1)*p_ie);
            rho_n = rho(1+e*p_ie);
            a = 0.5*(1/rho_n-1/rho_1);
            b = 0.5*(1/rho_n+1/rho_1);
            r = 1./(a*Q + b);
            J_1 = -a./(a*Q + b).^2;
            x = 1./r;
            x_i = 1./rho(1+(e-1)*p_ie:1+e*p_ie);
            [L,dLdx] = lagrangePolynomials(x,1:p_ie+1,p_ie+1,x_i);
            R = cell(1,2);
            R{1} = L;
            dxdr = -1./r.^2;
            dRdr = dLdx.*dxdr;
            J_2 = 1;
        else
            pts = controlPts(sctr,:);
            wgts = weights(element2(e,:),:);
            J_2 = prod(Xi_e(:,2)-Xi_e(:,1))/2;

            zeta = parent2ParametricSpace(Xi_e, Q);
            I = findKnotSpans(degree, zeta(1,:), knots);
            R = NURBSbasis(I, zeta, degree, knots, wgts);

            x = R{1}*pts;
            r = r_a./x;
            dxdzeta = R{2}*pts;
            J_1 = -r_a./x.^2.*dxdzeta; % = drdzeta
            dRdr = R{2}./J_1;
        end
        
        k1_e = zeros(n_en);
        k1_1_e = zeros(n_en);
        k1_2_e = zeros(n_en);
        k2_e = zeros(n_en);
        k3_e = zeros(n_en);
        k4_e = zeros(n_en);
        k5_e = zeros(n_en);
        for i = 1:numel(W)
            RR = R{1}(i,:).'*R{1}(i,:);
            RdRdr = R{1}(i,:).'*dRdr(i,:);
            fact = abs(J_1(i)) * J_2 * W(i);
            rUps = r(i)^2-Upsilon^2;
            kr = k*r(i);
            switch formulation
                case {'WBGC', 'WBGU'}
                    fact = fact * (r_a/r(i))^2;
            end
            switch formulation
                case {'PGC', 'BGC', 'WBGC'}
                    k1_e   = k1_e   + rUps*dRdr(i,:).'*dRdr(i,:) * fact;  
                    k1_1_e = k1_1_e + rUps*1i*k*(RdRdr - RdRdr.') * fact;
                    k1_2_e = k1_2_e - varrho3^2*RR                * fact;   
                case {'PGU', 'BGU', 'WBGU'}
                    fact = fact * exp(2*1i*kr);
                    k1_e = k1_e + (  rUps*dRdr(i,:).'*dRdr(i,:) ... 
                                   + rUps*1i*k*(RdRdr + RdRdr.') ...
                                   - (2*kr^2-varrho3^2)*RR) * fact;  
            end
            k2_e = k2_e +                  RR * fact;  
            k3_e = k3_e + varrho3^2      * RR * fact;  
            k4_e = k4_e + r(i)^2/rUps    * RR * fact;   
            k5_e = k5_e - Upsilon^2/rUps * RR * fact;  
        end
        spIdxRow(:,e) = kron(ones(1,n_en),sctr);
        spIdxCol(:,e) = kron(sctr,ones(1,n_en));
        
        K1values(:,e) = reshape(k1_e, sizeKe, 1);
        K1_1values(:,e) = reshape(k1_1_e, sizeKe, 1);
        K1_2values(:,e) = reshape(k1_2_e, sizeKe, 1);
        K2values(:,e) = reshape(k2_e, sizeKe, 1);
        K3values(:,e) = reshape(k3_e, sizeKe, 1);
        K4values(:,e) = reshape(k4_e, sizeKe, 1);
        K5values(:,e) = reshape(k5_e, sizeKe, 1);
    end
    Ntot = N-1 + noDofs_ie;
    K1 = shiftMatrix(Ntot-N,Ntot,K1);
    switch formulation
        case {'BGC', 'PGC', 'WBGC'}  
            K1_1 = shiftMatrix(Ntot-N,Ntot,K1_1);
            K1_2 = shiftMatrix(Ntot-N,Ntot,K1_2);
            
            K1_1 = K1_1 + sparse(spIdxRow,spIdxCol,K1_1values,Ntot,Ntot,numel(K1_1values));
            K1_2 = K1_2 + sparse(spIdxRow,spIdxCol,K1_2values,Ntot,Ntot,numel(K1_2values));
    end
    K2 = shiftMatrix(Ntot-N,Ntot,K2);
    K3 = shiftMatrix(Ntot-N,Ntot,K3);
    K4 = shiftMatrix(Ntot-N,Ntot,K4);
    K5 = shiftMatrix(Ntot-N,Ntot,K5);
    
    K1 = K1 + sparse(spIdxRow,spIdxCol,K1values,Ntot,Ntot,numel(K1values));
    K2 = K2 + sparse(spIdxRow,spIdxCol,K2values,Ntot,Ntot,numel(K2values));
    K3 = K3 + sparse(spIdxRow,spIdxCol,K3values,Ntot,Ntot,numel(K3values));
    K4 = K4 + sparse(spIdxRow,spIdxCol,K4values,Ntot,Ntot,numel(K4values));
    K5 = K5 + sparse(spIdxRow,spIdxCol,K5values,Ntot,Ntot,numel(K5values));
    
    noDofs_new = noDofs_new + noSurfDofs*(noDofs_ie-1);
else
    Ntot = N;
end
switch formulation
    case {'PGU', 'BGU', 'WBGU'}
        varrho2 = k*r_a;
        K1 = K1*exp(-2*1i*varrho2);
        K2 = K2*exp(-2*1i*varrho2);
        K3 = K3*exp(-2*1i*varrho2);
        K4 = K4*exp(-2*1i*varrho2);
        K5 = K5*exp(-2*1i*varrho2);
end
if task.rom.useROM
    A_2 =  kron(K1_2,A1values) ...
         + kron(K3,A3values);
    A_1 =  kron(K1_1,A1values);
    A =    kron(K1,A1values) ...
         + kron(K2,A2values) ...
         + kron(K4,A4values) ...
         + kron(K5,A5values);
else
    switch formulation
        case {'PGC', 'BGC', 'WBGC'}
            K1 = K1 + K1_1*omega + K1_2*omega^2;
    end
    A =   kron(K1,A1values) ...
        + kron(K2,A2values) ...
        + kron(K3,A3values) ...
        + kron(K4,A4values) ...
        + kron(K5,A5values);
end

dofsToRemove = setdiff(1:noSurfDofs,unique(spIdx(:,1)));
noDofsToRemove = numel(dofsToRemove);
newDofsToRemove = zeros(1,noDofsToRemove*Ntot);
i = 1;
for n = 1:Ntot
    newDofsToRemove(i:i+noDofsToRemove-1) = dofsToRemove+(n-1)*noSurfDofs;
    i = i + noDofsToRemove;
end

[i,j,Avalues] = find(A);    
if noDofs > 0
    [i,j,newDofsToRemove] = modifyIndices(i,j,noSurfDofs,noDofs,nodes,newDofsToRemove);    
    task.varCol{1}.Ainf = sparse(i,j,Avalues,noDofs_new,noDofs_new);
    
    if task.rom.useROM
        [i,j,Avalues] = find(A_1);
        [i,j] = modifyIndices(i,j,noSurfDofs,noDofs,nodes);
        task.varCol{1}.Ainf1 = sparse(i,j,Avalues,noDofs_new,noDofs_new);
        
        [i,j,Avalues] = find(A_2);
        [i,j] = modifyIndices(i,j,noSurfDofs,noDofs,nodes);
        task.varCol{1}.Ainf2 = sparse(i,j,Avalues,noDofs_new,noDofs_new);
    end
else
    task.varCol{1}.Ainf = sparse(i,j,Avalues,noDofs_new,noDofs_new);
end
dofsToRemove_old = task.varCol{1}.dofsToRemove;
task.varCol{1}.dofsToRemove = sort(unique([dofsToRemove_old newDofsToRemove]));
task.varCol{1}.dofsToRemove_old = dofsToRemove_old;
task.varCol{1}.noDofs_new = noDofs_new;

function A = shiftMatrix(shift,noDofs,A)
[i,j,Kvalues] = find(A);
i = i+shift;
j = j+shift;
A = sparse(i,j,Kvalues,noDofs,noDofs);

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



