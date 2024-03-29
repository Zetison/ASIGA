function task = infElementsNonSepGeom(task)

degree = task.varCol{1}.degree(1:2); % assume degree is equal in all patches

elRange = task.varCol{1}.elRange;
weights = task.varCol{1}.weights;
controlPts = task.varCol{1}.controlPts;
knotVecs = task.varCol{1}.knotVecs;

extraGP = task.misc.extraGP;

N = task.iem.N;
Ntot = task.iem.N;
formulation = task.misc.formulation;

noDofs = task.varCol{1}.noDofs;

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

%Use chebychev polynomials
chimin = task.varCol{1}.chimin;
chimax = task.varCol{1}.chimax;
x_0 = task.iem.x_0;
A_2 = task.iem.A_2;
[D,Dt] = generateCoeffMatrix(task);

k = task.misc.omega/task.varCol{1}.c_f;
Upsilon = task.iem.Upsilon;

useApproxRadialIntegrals = 2; % 0 is most accurate, 1 is a spline approximation (not reccomended here as Chebychev is here implemented more robust), 2 is Chebychev (reccomended)
B1splines = {};
B2splines = {};
B1cheb = {};
B2cheb = {};

if useApproxRadialIntegrals == 1
    p_splineApprox = 2;
    B1splines = cell(2*N+4,1);
    parfor i = 1:2*N+4
        B1splines{i} = splineInterpolation(@(chi)radialIntegral2(i, chi, k, Upsilon, formulation, 1),...
                                           chimin,chimax,max(ceil((chimax-chimin)*k/3*200),p_splineApprox+1),p_splineApprox);
    end
    B2splines = cell(2*N+3,1);
    parfor i = 1:2*N+3
        B2splines{i} = splineInterpolation(@(chi)radialIntegral2(i, chi, k, Upsilon, formulation, 2),...
                                           chimin,chimax,max(ceil((chimax-chimin)*k/3*200),p_splineApprox+1),p_splineApprox);
    end
elseif useApproxRadialIntegrals == 2
    B1cheb = cell(2*N+4,1);
    parfor i = 1:2*N+4
        B1cheb{i} = adaptiveChebychevInterp(@(x)radialIntegral2(i, chimin+0.5*(chimax-chimin)*(1+x), k, Upsilon, ...
                                                formulation, 1), 1e-6);
    end
    B2cheb = cell(2*N+3,1);
    parfor i = 1:2*N+3
        B2cheb{i} = adaptiveChebychevInterp(@(x)radialIntegral2(i, chimin+0.5*(chimax-chimin)*(1+x), k, Upsilon, ...
                                                formulation, 2), 1e-6);
    end
end

%% Calculate contribution from infinite elements
spIdxRow = zeros((N*n_en)^2,noElems);
spIdxCol = zeros((N*n_en)^2,noElems);
Avalues = zeros((N*n_en)^2,noElems);
[Q, W] = gaussTensorQuad(degree+1+extraGP(1:2));


progressBars = task.misc.progressBars;
nProgressStepSize = ceil(noElems/1000);
if progressBars
    try
        ppm = ParforProgMon('Building IE matrix: ', noElems, nProgressStepSize);
    catch
        progressBars = false;
        ppm = NaN;
    end
else
    ppm = NaN;
end

max_r_a_recorded = -inf;
min_r_a_recorded = inf;
% totArea = 0;
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
    wgts = weights(element2(e,:)); % New  
            
    Rs = NURBSbasis(I, xi, degree, knots, wgts);
    dXdxi = Rs{2}*pts;
    dXdeta = Rs{3}*pts;

    X = Rs{1}*pts;
    Xt = (X-x_0)*A_2.';

    [r_a_arr, theta_arr, ~, c1_arr, c2_arr] = evaluateProlateCoords(Xt,Upsilon); 
    if useApproxRadialIntegrals == 1
        B1 = zeros(size(W,1),2*N+4);
        xi_pt = (r_a_arr-chimin)/(chimax-chimin);
        for i = 1:2*N+4
            for gp = 1:size(W,1)
                B1(gp,i) = evaluateNURBS(B1splines{i}{1}, xi_pt(gp));
            end
        end
        B2 = zeros(size(W,1),2*N+3);
        for i = 1:2*N+3
            for gp = 1:size(W,1)
                B2(gp,i) = evaluateNURBS(B2splines{i}{1}, xi_pt(gp));
            end
        end
    elseif useApproxRadialIntegrals == 2
        B1 = zeros(size(W,1),2*N+4);
        x_cheb = (2*r_a_arr-chimax-chimin)/(chimax-chimin);

        for i = 1:2*N+4
            for j = 1:length(B1cheb{i})
                if j == 1 || j == length(B1cheb{i})
                    B1(:,i) = B1(:,i) + 0.5*B1cheb{i}(j)*chebyshevTilde(j-1,x_cheb);
                else
                    B1(:,i) = B1(:,i) + B1cheb{i}(j)*chebyshevTilde(j-1,x_cheb);
                end
            end
        end
        B2 = zeros(size(W,1),2*N+3);
        for i = 1:2*N+3
            for j = 1:length(B2cheb{i})
                if j == 1 || j == length(B2cheb{i})
                    B2(:,i) = B2(:,i) + 0.5*B2cheb{i}(j)*chebyshevTilde(j-1,x_cheb);
                else
                    B2(:,i) = B2(:,i) + B2cheb{i}(j)*chebyshevTilde(j-1,x_cheb);
                end
            end
        end      
    else
        B1 = zeros(size(W,1),2*N+4);
        for i = 1:2*N+4
            B1(:,i) = radialIntegral2(i, r_a_arr, k, Upsilon, formulation, 1);
        end
        B2 = zeros(size(W,1),2*N+3);
        for i = 1:2*N+3
            B2(:,i) = radialIntegral2(i, r_a_arr, k, Upsilon, formulation, 2);
        end      
    end
    
    temp = zeros(n_en,n_en,N,N);
    for gp = 1:size(W,1)
        R = Rs{1}(gp,:);
        r_a = r_a_arr(gp);
        theta = theta_arr(gp);
        c1 = c1_arr(gp);
        c2 = c2_arr(gp);
        
        wt = W(gp);
        
        J = [dXdxi(gp,:).',dXdeta(gp,:).'];
                
        if r_a > chimax || r_a < chimin
            error('chi is not in range [chimin chimax]')
        end
%         if r_a < min_r_a_recorded
%             min_r_a_recorded = r_a;
%         end
%         if r_a > max_r_a_recorded
%             max_r_a_recorded = r_a;
%         end
        DPDX = dPdX(Xt(gp,:),Upsilon,r_a,c1,c2)*A_2;
        
        J3 = DPDX(2:3,:)*J;
        J_3 = abs(det(J3));
        
        
        [dr_adtheta, dr_adphi] = dchidP(DPDX(1,:),J,J3);
        
        dRdP = J3'\[Rs{2}(gp,:); Rs{3}(gp,:)];
        
        dRdtheta = dRdP(1,:);
        dRdphi = dRdP(2,:);

        RR = R'*R;
        RtRt = dRdtheta'*dRdtheta;  
        RpRp = dRdphi'*dRdphi;  
        
        RtR = dRdtheta'*R;  
        RRt = R'*dRdtheta; 
        RpR = dRdphi'*R; 
        RRp = R'*dRdphi;

        varrho1 = Upsilon/r_a;
        varrho2 = k*r_a;
        varrho3 = k*Upsilon;
        for nt = 1:N
            for mt = 1:N                
                switch formulation
                    case 'PGU'
                        temp2 = RR*(-2*varrho2^2*B1(gp,nt+mt) - 1i*varrho2*(nt+mt+2)*B1(gp,nt+mt+1) +((nt+2)*mt + varrho3^2*(1+cos(theta)^2))*B1(gp,nt+mt+2) ...
                                    +1i*varrho1*varrho3*(nt+mt+2)*B1(gp,nt+mt+3)-varrho1^2*(nt+2)*mt*B1(gp,nt+mt+4) ...
                                    +(-varrho2^2-1i*varrho2*(nt+mt+2)+(nt+2)*mt) ...
                                        *1/r_a^2*(dr_adtheta^2*B1(gp,nt+mt+2)+dr_adphi^2/sin(theta)^2*(B2(gp,nt+mt+1)-varrho1^2*cos(theta)^2*B2(gp,nt+mt+3)))) ...
                                +(RtRt + (mt-1i*varrho2)/r_a*dr_adtheta*RtR + (nt+2-1i*varrho2)/r_a*dr_adtheta*RRt)*B1(gp,nt+mt+2) ...
                                +(RpRp + (mt-1i*varrho2)/r_a*dr_adphi*RpR + (nt+2-1i*varrho2)/r_a*dr_adphi*RRp) ...
                                    *1/sin(theta)^2*(B2(gp,nt+mt+1) - varrho1^2*cos(theta)^2*B2(gp,nt+mt+3));
                    case 'PGC'
                        temp2 = RR*(- 1i*varrho2*(nt-mt+2)*B1(gp,nt+mt+1) + ((nt+2)*mt + varrho3^2*(-1+cos(theta)^2))*B1(gp,nt+mt+2) ...
                                    +1i*varrho1*varrho3*(nt-mt+2)*B1(gp,nt+mt+3)-varrho1^2*(nt+2)*mt*B1(gp,nt+mt+4) ...
                                    +(varrho2^2-1i*varrho2*(nt-mt+2)+(nt+2)*mt) ...
                                        *1/r_a^2*(dr_adtheta^2*B1(gp,nt+mt+2)+dr_adphi^2/sin(theta)^2*(B2(gp,nt+mt+1)-varrho1^2*cos(theta)^2*B2(gp,nt+mt+3)))) ...
                                +(RtRt + (mt-1i*varrho2)/r_a*dr_adtheta*RtR + (nt+2+1i*varrho2)/r_a*dr_adtheta*RRt)*B1(gp,nt+mt+2) ...
                                +(RpRp + (mt-1i*varrho2)/r_a*dr_adphi*RpR + (nt+2+1i*varrho2)/r_a*dr_adphi*RRp) ...
                                    *1/sin(theta)^2*(B2(gp,nt+mt+1) - varrho1^2*cos(theta)^2*B2(gp,nt+mt+3));
                    case 'BGU'
                        if mt+nt == 2
                            temp2 = RR*(-2*1i*varrho2*B1(gp,1) +(1 + varrho3^2*(1+cos(theta)^2))*B1(gp,2) ...
                                        +2*1i*varrho1*varrho3*B1(gp,3)-varrho1^2*B1(gp,4) - 1i*varrho2*exp(2*1i*varrho2) ...
                                        +(-varrho2^2-2*1i*varrho2+1) ...
                                            *1/r_a^2*(dr_adtheta^2*B1(gp,2)+dr_adphi^2/sin(theta)^2*(B2(gp,1)-varrho1^2*cos(theta)^2*B2(gp,3)))) ...
                                    +(RtRt + (1-1i*varrho2)/r_a*dr_adtheta*RtR + (1-1i*varrho2)/r_a*dr_adtheta*RRt)*B1(gp,2) ...
                                    +(RpRp + (1-1i*varrho2)/r_a*dr_adphi*RpR + (1-1i*varrho2)/r_a*dr_adphi*RRp) ...
                                        *1/sin(theta)^2*(B2(gp,1) - varrho1^2*cos(theta)^2*B2(gp,3));
                        else
                            temp2 = RR*(-2*varrho2^2*B1(gp,nt+mt-2) - 1i*varrho2*(nt+mt)*B1(gp,nt+mt-1) +(nt*mt + varrho3^2*(1+cos(theta)^2))*B1(gp,nt+mt) ...
                                        +1i*varrho1*varrho3*(nt+mt)*B1(gp,nt+mt+1)-varrho1^2*nt*mt*B1(gp,nt+mt+2) ...
                                        +(-varrho2^2-1i*varrho2*(nt+mt)+nt*mt) ...
                                            *1/r_a^2*(dr_adtheta^2*B1(gp,nt+mt)+dr_adphi^2/sin(theta)^2*(B2(gp,nt+mt-1)-varrho1^2*cos(theta)^2*B2(gp,nt+mt+1)))) ...
                                    +(RtRt + (mt-1i*varrho2)/r_a*dr_adtheta*RtR + (nt-1i*varrho2)/r_a*dr_adtheta*RRt)*B1(gp,nt+mt) ...
                                    +(RpRp + (mt-1i*varrho2)/r_a*dr_adphi*RpR + (nt-1i*varrho2)/r_a*dr_adphi*RRp) ...
                                        *1/sin(theta)^2*(B2(gp,nt+mt-1) - varrho1^2*cos(theta)^2*B2(gp,nt+mt+1));
                        end
                    case 'BGC'
                        if mt+nt == 2
                            temp2 = RR*((1+varrho3^2*(-1+cos(theta)^2))*B1(gp,2) ...
                                        -varrho1^2*B1(gp,4) - 1i*varrho2...
                                        +(varrho2^2+1) ...
                                            *1/r_a^2*(dr_adtheta^2*B1(gp,2)+dr_adphi^2/sin(theta)^2*(B2(gp,1)-varrho1^2*cos(theta)^2*B2(gp,3)))) ...
                                    +(RtRt + (1-1i*varrho2)/r_a*dr_adtheta*RtR+(1+1i*varrho2)/r_a*dr_adtheta*RRt)*B1(gp,2) ...
                                    +(RpRp + (1-1i*varrho2)/r_a*dr_adphi*RpR+(1+1i*varrho2)/r_a*dr_adphi*RRp) ...
                                        *1/sin(theta)^2*(B2(gp,1) - varrho1^2*cos(theta)^2*B2(gp,3));
                        else
                            temp2 = RR*(-1i*varrho2*(nt-mt)*B1(gp,nt+mt-1) + (nt*mt+varrho3^2*(-1+cos(theta)^2))*B1(gp,nt+mt) ...
                                        +1i*varrho1*varrho3*(nt-mt)*B1(gp,nt+mt+1)-varrho1^2*nt*mt*B1(gp,nt+mt+2) ...
                                        +(varrho2^2-1i*varrho2*(nt-mt)+nt*mt) ...
                                            *1/r_a^2*(dr_adtheta^2*B1(gp,nt+mt)+dr_adphi^2/sin(theta)^2*(B2(gp,nt+mt-1)-varrho1^2*cos(theta)^2*B2(gp,nt+mt+1)))) ...
                                    +(RtRt + (mt-1i*varrho2)/r_a*dr_adtheta*RtR+(nt+1i*varrho2)/r_a*dr_adtheta*RRt)*B1(gp,nt+mt) ...
                                    +(RpRp + (mt-1i*varrho2)/r_a*dr_adphi*RpR+(nt+1i*varrho2)/r_a*dr_adphi*RRp) ...
                                        *1/sin(theta)^2*(B2(gp,nt+mt-1) - varrho1^2*cos(theta)^2*B2(gp,nt+mt+1));
                        end
                end
                for n = 1:N
                    for m = 1:N
                        if Dt(n,nt)*D(m,mt) ~= 0
                            switch formulation
                                case {'PGU', 'BGU'}
                                    temp(:,:,n,m) = temp(:,:,n,m) + temp2*Dt(n,nt)*D(m,mt)*r_a*exp(-2*1i*varrho2)*sin(theta)*J_3*J_2*wt;
                                case {'PGC', 'BGC'}
                                    temp(:,:,n,m) = temp(:,:,n,m) + temp2*Dt(n,nt)*D(m,mt)*r_a*sin(theta)*J_3*J_2*wt;
                            end
                        end
                    end
                end
            end
        end
%         totArea = totArea + J_3*J_2*wt;
    end  
    A_inf_values_temp = zeros((N*n_en)^2,1);
    spIdxRow_temp = zeros((N*n_en)^2,1);
    spIdxCol_temp = zeros((N*n_en)^2,1);
    counter = 1;
    for n = 1:N
        for m = 1:N
            indices = counter:(counter+n_en^2-1);
            A_inf_values_temp(indices) = reshape(temp(:,:,n,m),n_en^2,1);
            spIdxRow_temp(indices) = copyVector(sctrLocal+(noSurfDofs*(n-1)),n_en,1);
            spIdxCol_temp(indices) = copyVector(sctrLocal+(noSurfDofs*(m-1)),n_en,2);
            counter = counter + n_en^2;
        end
    end
    Avalues(:,e) = A_inf_values_temp;
    
    spIdxRow(:,e) = spIdxRow_temp;
    spIdxCol(:,e) = spIdxCol_temp;
end
% errorInTotArea = totArea-2*pi^2;

% max_r_a_recorded
% min_r_a_recorded

spIdxRow = reshape(spIdxRow,numel(spIdxRow),1);
spIdxCol = reshape(spIdxCol,numel(spIdxCol),1);
Avalues = reshape(Avalues,numel(Avalues),1);

[spIdx,~,IuniqueIdx] = unique([spIdxRow, spIdxCol],'rows');
Avalues = accumarray(IuniqueIdx,Avalues);

dofsToRemove = setdiff(1:noSurfDofs,unique(spIdx(:,1)));
noDofsToRemove = numel(dofsToRemove);
newDofsToRemove = zeros(1,noDofsToRemove*Ntot);
i = 1;
for n = 1:Ntot
    newDofsToRemove(i:i+noDofsToRemove-1) = dofsToRemove+(n-1)*noSurfDofs;
    i = i + noDofsToRemove;
end
ii = spIdx(:,1);
jj = spIdx(:,2);
if noDofs > 0
    [ii,jj,newDofsToRemove] = modifyIndices(ii,jj,noSurfDofs,noDofs,nodes,newDofsToRemove);    
    task.varCol{1}.Ainf = sparse(ii,jj,Avalues,noDofs_new,noDofs_new);
else
    task.varCol{1}.Ainf = sparse(ii,jj,Avalues,noDofs_new,noDofs_new,numel(IuniqueIdx));
end
dofsToRemove_old = task.varCol{1}.dofsToRemove;
task.varCol{1}.dofsToRemove = sort(unique([dofsToRemove_old newDofsToRemove]));
task.varCol{1}.dofsToRemove_old = dofsToRemove_old;
task.varCol{1}.noDofs_new = noDofs_new;


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



