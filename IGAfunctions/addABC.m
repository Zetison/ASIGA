function task = addABC(task)

degree = task.varCol{1}.degree(1:2); % assume degree is equal in all patches

N = task.abc.N;

k = task.varCol{1}.k;
Upsilon = task.iem.Upsilon;
r_a = task.misc.r_a;

formulation = task.misc.formulation;
A_2 = task.iem.A_2;
x_0 = task.iem.x_0;
noDofs = task.varCol{1}.noDofs;

weights = task.varCol{1}.weights;
controlPts = task.varCol{1}.controlPts;

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

%% Calculate contribution from Gamma_a
spIdxRow = zeros(n_en^2,noElems);
spIdxCol = zeros(n_en^2,noElems);
Avalues = zeros(n_en^2,noElems);
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
    
    A_IJ = zeros(n_en);
            
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
        wt = W(gp);
        theta = theta_arr(gp);
        RR = R{1}(gp,:)'*R{1}(gp,:);

        DPDX = dPdX(Xt(gp,:),Upsilon,r_a,d1(gp),d2(gp))*A_2;
        J = [a11(gp) a12(gp);
             a21(gp) a22(gp);
             a31(gp) a32(gp)];
             
        J3 = DPDX(2:3,:)*J;
        J_3 = abs(det(J3));

        dRdP = J3'\[R{2}(gp,:); R{3}(gp,:)];
        dRdtheta = dRdP(1,:);
        dRdphi = dRdP(2,:);

        switch formulation
            case 'HH'
                switch N
                    case 1
                        A_IJ = A_IJ - (1i*k-1/r_a)*RR*sqrt(r_a^2-Upsilon^2*cos(theta)^2)*sqrt(r_a^2-Upsilon^2)*sin(theta)*J_3*J_2*wt;   
                    case 2
                        error('Not properly implemented')
                        A_IJ = A_IJ - ((1i*k-1/r_a)*RR ...
                                       - 1/(2*r_a^3)/(1-1i*k*r_a)*(dRdtheta'*dRdtheta + dRdphi'*dRdphi/sin(theta)^2)) ...
                                      *sqrt(r_a^2-Upsilon^2*cos(theta)^2)*sqrt(r_a^2-Upsilon^2)*sin(theta)*J_3*J_2*wt; 
    %                         A_IJ = A_IJ + ((1i*k-1/r_a)*RR ...
    %                                        - 1/(2*r_a)/(1-1i*k*r_a)*(R'*dR2dtheta2 + cot(theta)*R'*dRdtheta + R'*dR2dphi2/sin(theta)^2)) ...
    %                                       *sqrt(r_a^2-Upsilon^2*cos(theta)^2)*sqrt(r_a^2-Upsilon^2)*sin(theta)*J_3*J_2*wt;  
                end
            case 'BGT'
                switch N
                    case 1
                        A_IJ = A_IJ - (1i*k-1/r_a)*RR*sqrt(r_a^2-Upsilon^2*cos(theta)^2)*sqrt(r_a^2-Upsilon^2)*sin(theta)*J_3*J_2*wt;  
                    case 2
                        error('Not properly implemented')
                        A_IJ = A_IJ - ((1i*k-1/r_a)*RR ...
                                       - 1/(2*r_a^3)/(1-1i*k*r_a)*(dRdtheta'*dRdtheta + dRdphi'*dRdphi/sin(theta)^2)) ...
                                      *sqrt(r_a^2-Upsilon^2*cos(theta)^2)*sqrt(r_a^2-Upsilon^2)*sin(theta)*J_3*J_2*wt; 
    %                         A_IJ = A_IJ + ((1i*k-1/r_a)*RR ...
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

task.varCol{1}.Ainf = sparse(spIdx(:,1),spIdx(:,2),Avalues,noDofs,noDofs,numel(IuniqueIdx));
