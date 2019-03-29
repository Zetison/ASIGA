function [A_inf, newDofsToRemove] = addInfElements2D(varCol, k, f, r_a)

elRangeXi = varCol.elRangeXi;
Xi = varCol.nurbs.knots{1};

p = varCol.nurbs.degree(1);

n = varCol.nurbs.number(1);
m = varCol.nurbs.number(2);

gluedNodes = varCol.gluedNodes;
N = varCol.N;
infiniteElementFormulation = varCol.infiniteElementFormulation;
rm = varCol.rm;

noDofs = varCol.noDofs;

weights = varCol.weights;
controlPts = varCol.controlPts;
nurbs = varCol.nurbs;

eta1Nodes = zeros(1,n);
counter = 1;
for i = 1:n
    eta1Nodes(counter) = n*(m-1) + i;
    counter = counter + 1;
end

[XiMesh, indexXi, noElemsXi] = generateIGA1DMesh(Xi, p, n);

% Glue nodes in 1D mesh
for i = 1:length(gluedNodes)
    parentIdx = gluedNodes{i}(1);
    for j = 2:length(gluedNodes{i})
        indices = (XiMesh == gluedNodes{i}(j));
        XiMesh(indices) = parentIdx;
    end
end

%% Find coefficients for radial shape functions in infinite elements
S1 = zeros(N);
S2 = zeros(N);
for j = 1:N
    for jm = 1:N
        switch infiniteElementFormulation
            case {'PGC', 'PGU'}
                n = j+2;
            case {'BGC', 'BGU'}
                n = j;
        end
        m = j;
        S1(j,jm) = rm(jm)^(-n+1/2);
        S2(j,jm) = rm(jm)^(-m+1/2);
    end
end
H1_temp = inv(S1);
H2_temp = inv(S2);


H1 = zeros(N);
H2 = zeros(N);
for i = 1:N
    for j = 1:N
        H1(i,j) = H1_temp(i,j);
        H2(i,j) = H2_temp(i,j);
    end
end

%% Evaluate analytic integrals in ``radial'' direction. 
% Note that the last two integrals I1(end) and I1(end-1),
% will be redundant for the cases 'BGC' and 'BGU'
I1 = zeros(2*N+4,1);
for n = 1:2*N+4
    I1(n) = radialIntegral(n-1, r_a, k, f, 2, infiniteElementFormulation, 1);
end
I2 = radialIntegral(NaN, r_a, k, f, 2, infiniteElementFormulation, 2);

%% Calculate contribution from infinite elements
n_en = p+1;

spIdxRow = zeros((N*n_en)^2,noElemsXi);
spIdxCol = zeros((N*n_en)^2,noElemsXi);
A_inf_values = zeros((N*n_en)^2,noElemsXi);

[W1D,Q1D] = gaussianQuadNURBS(p+1); 
% for e = 1:noElemsXi
parfor e = 1:noElemsXi
    idXi = indexXi(e,1); 

    Xi_e = elRangeXi(idXi,:);

    J_2 = 0.5*(Xi_e(2)-Xi_e(1));

    sctrXi = eta1Nodes(XiMesh(e,:));
    
    pts = controlPts(sctrXi,:);
    
    J1_AB = zeros(n_en);
    J2_AB = zeros(n_en);
    J3_AB = zeros(n_en);
    
    for gp = 1:size(W1D,1)
        pt = Q1D(gp,:);
        wt = W1D(gp);

        xi  = parent2ParametricSpace(Xi_e, pt(1));
        
        [R_fun, dRdxi] = NURBS1DBasis(xi, p, Xi, weights(eta1Nodes));
        J = pts'*dRdxi';
        
        v = evaluateNURBS(nurbs,[xi 1]);
        x = v(1);
        y = v(2);
        
        [radius, theta] = evaluateEllipticCoords(x,y,f);
        
        DPDX = dPdX_elliptic(f,radius,theta);
        
        J3 = DPDX(2,:)*J;
        J_3 = abs(J3);
        
        dRdtheta = dRdxi/J3;
        
        J1_AB = J1_AB  + R_fun'*R_fun                                  *J_3*J_2*wt;  
        J2_AB = J2_AB  + dRdtheta'*dRdtheta                            *J_3*J_2*wt;  
        J3_AB = J3_AB  + R_fun'*R_fun*cos(theta)^2                     *J_3*J_2*wt;   
%         keyboard
    end   
    
    A_inf_values_temp = zeros((N*n_en)^2,1);
    spIdxRow_temp = zeros((N*n_en)^2,1);
    spIdxCol_temp = zeros((N*n_en)^2,1);
    for j = 1:N
        for i = 1:N
            switch infiniteElementFormulation
                case {'PGC', 'PGU'}
                    n = i + 2;
                case {'BGC', 'BGU'}
                    n = i;
            end
            m = j;
            temp = zeros(n_en);
            for im = 1:N
                for jm = 1:N
                    switch infiniteElementFormulation
                        case {'PGC', 'PGU'}
                            nm = im + 2;
                        case {'BGC', 'BGU'}
                            nm = im;
                    end
                    mm = jm;
                    switch infiniteElementFormulation
                        case {'PGU', 'BGU'}
                            if mm == 1 && nm == 1
                                blf  = J1_AB*((1/4+k^2*f^2)*I1(2) - 1i*k*(I1(1) - f^2*I1(3)) - f^2/4*I1(4)) ...
                                       + J2_AB*I1(2) + J3_AB*k^2*f^2*I1(2) - 1i*k*J1_AB*(r_a*exp(2*1i*k*r_a)/sqrt(r_a^2-f^2)-f^2*I2);
                            else
                                blf  = J1_AB*(-2*k^2*I1(nm+mm-2) + k^2*f^2*I1(nm+mm) + 1i*k*(1-nm-mm)*(I1(nm+mm-1)-f^2*I1(nm+mm+1))...
                                       +(2*nm-1)*(2*mm-1)/4*(I1(nm+mm)-f^2*I1(nm+mm+2))) + J2_AB*I1(nm+mm) + J3_AB*k^2*f^2*I1(nm+mm);
                            end
                        case {'PGC', 'BGC'}
                            if mm == 1 && nm == 1
                                slf  = J1_AB*((1/4-k^2*f^2)*I1(2) - f^2/4*I1(4) - 1i*k) ...
                                       + J2_AB*I1(2) + J3_AB*k^2*f^2*I1(2);
                            else
                                slf  = J1_AB*(-k^2*f^2*I1(nm+mm) + 1i*k*(mm-nm)*(I1(nm+mm-1) - f^2*I1(nm+mm+1)) ...
                                        +(2*nm-1)*(2*mm-1)/4*(I1(nm+mm) - f^2*I1(nm+mm+2))) + J2_AB*I1(nm+mm) + J3_AB*k^2*f^2*I1(nm+mm);
                            end
                    end
                    switch infiniteElementFormulation
                        case {'PGU', 'BGU'}
                            temp = temp + H1(i,im)*exp(-1i*k*rm(i))*H2(j,jm)*exp(-1i*k*rm(j))*blf;
                        case {'PGC', 'BGC'}
                            temp = temp + H1(i,im)*exp(-1i*k*rm(i))*H2(j,jm)*exp(1i*k*rm(j))*slf;
                    end
                end
            end
            indices = (1:n_en^2)+(n_en^2*N*(j-1) + n_en^2*(i-1));
            A_inf_values_temp(indices) = reshape(temp,n_en^2,1);
            spIdxRow_temp(indices) = copyVector(sctrXi+(noDofs*(i-1)),n_en,1);
            spIdxCol_temp(indices) = copyVector(sctrXi+(noDofs*(j-1)),n_en,2);
        end
    end
    A_inf_values(:,e) = A_inf_values_temp;
    spIdxRow(:,e) = spIdxRow_temp;
    spIdxCol(:,e) = spIdxCol_temp;
end

noDofs_new = noDofs + noDofs*(N-1);

A_inf = sparse(spIdxRow,spIdxCol,A_inf_values);


if min(size(A_inf)) < noDofs_new
    A_inf(noDofs+(N-1)*noDofs, noDofs+(N-1)*noDofs) = 0;
end

newDofsToRemove = setdiff((noDofs+1):noDofs_new,unique(spIdxRow));


