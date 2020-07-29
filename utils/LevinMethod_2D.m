function I = LevinMethod_2D(f,g,dgdx,dgdy,omega)
useGalerkin = 1;
spline = getRectangleData(2, 2);
degreeElev = 2;
newXiKnots = 0;
newEtaKnots = newXiKnots;

spline = elevateNURBSdegree(spline,degreeElev*[1, 1]);
spline = insertKnotsInNURBS(spline,{insertUniform2(spline.knots{1}, newXiKnots)
                                    insertUniform2(spline.knots{2}, newEtaKnots)});

varCol.dimension = 1;    
varCol = convertNURBS(spline, varCol);
varCol = generateIGA2DMesh_new(varCol);

p_xi = spline.degree(1);
p_eta = spline.degree(2);
noElems = varCol.noElems;
controlPts = varCol.controlPts;
noDofs = spline.number(1)*spline.number(2);

element = varCol.element;
weights = varCol.weights;
Xi = spline.knots{1};
Eta = spline.knots{2};

if useGalerkin
    elRangeXi = varCol.elRangeXi;
    elRangeEta = varCol.elRangeEta;
    index = varCol.index;
    [W2D,Q2D] = gaussianQuadNURBS(p_xi+1,p_eta+1);
    n_en = (p_xi+1)*(p_eta+1);
    Avalues = zeros(n_en^2,noElems); 
    spIdxRow = zeros(n_en^2,noElems);
    spIdxCol = zeros(n_en^2,noElems);
    F_indices = zeros(n_en,noElems); 
    Fvalues = zeros(n_en,noElems); 
    parfor e = 1:noElems
        idXi = index(e,1);
        idEta = index(e,2);

        Xi_e = elRangeXi(idXi,:);
        Eta_e = elRangeEta(idEta,:);

        sctr = element(e,:);

        pts = controlPts(sctr,:);
        a_e = zeros(n_en);
        f_e = zeros(n_en,1);
        
        J_2 = 0.25*(Xi_e(2)-Xi_e(1))*(Eta_e(2)-Eta_e(1));
        for gp = 1:size(W2D,1)
            pt = Q2D(gp,:);
            wt = W2D(gp);

            xi  = parent2ParametricSpace(Xi_e, pt(1));
            eta  = parent2ParametricSpace(Eta_e, pt(2));

            [R, dRdxi, dRdeta] = NURBS2DBasis2(xi, eta, p_xi, p_eta, Xi, Eta, weights);
            J = pts'*[dRdxi' dRdeta'];
            dRdX = J'\[dRdxi; dRdeta];
            J_1 = det(J);
            x = R*pts;
            a_e = a_e + R'*(dRdX(1,:) - dRdX(2,:) + 1i*omega*(dgdx(x(1),x(2)) - dgdy(x(1),x(2)))*R)* abs(J_1) * J_2 * wt;  
            f_e = f_e + R'*f(x(1),x(2)) * abs(J_1) * J_2 * wt;
        end
        spIdxRow(:,e) = copyVector(sctr,n_en,1);
        spIdxCol(:,e) = copyVector(sctr,n_en,2);
        Avalues(:,e) = reshape(a_e, n_en^2, 1);
        F_indices(:,e) = sctr;
        Fvalues(:,e) = f_e;
    end
    A = sparse(spIdxRow,spIdxCol,Avalues);
    FF = vectorAssembly(Fvalues,F_indices,noDofs);
    if 0
        A(1,:) = [];
        A(:,1) = [];
        FF(1) = [];
        U = A\FF;
        U = [0; U];
    else
        U = A\FF;
    end
    cond(full(A))
else
    uniqueXi = unique(Xi);
    uniqueEta = unique(Eta);

    grev_xi = aveknt(Xi, p_xi+1);
    grev_eta = aveknt(Eta, p_eta+1);
    n = length(grev_xi);
    cp = zeros(n^2,2);
    counter = 1;
    for i = 1:n
        for j = 1:n
            cp(counter,:) = [grev_xi(i), grev_eta(j)];
            counter = counter + 1;
        end
    end
    noElementsXi = length(uniqueXi) - 1;
    noElementsEta = length(uniqueEta) - 1;

    A = zeros(noDofs);
    x_arr = zeros(noDofs,2);
    for i = 1:length(cp)
        xi = cp(i,1);
        eta = cp(i,2);
        xi_idx = findKnotSpan(noElementsXi, 0, xi, uniqueXi);
        eta_idx = findKnotSpan(noElementsEta, 0, eta, uniqueEta);
        e = xi_idx + noElementsXi*(eta_idx-1);
        sctr = element(e,:);

        [R, dRdxi, dRdeta] = NURBS2DBasis(xi, eta, p_xi, p_eta, Xi, Eta, weights);
        pts = controlPts(sctr,:);
        J_1 = pts'*[dRdxi' dRdeta'];
        dRdX = J_1'\[dRdxi; dRdeta];
        x = R*pts;
        x_arr(i,:) = x;
        A(i,sctr) = dRdX(1,:) - dRdX(2,:) + 1i*omega*(dgdx(x(1),x(2)) - dgdy(x(1),x(2)))*R;
    end

    FF = f(x_arr(:,1), x_arr(:,2));
    U = A\FF;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if ~useGalerkin
%     cp_x = x_arr(:,1);
%     cp_y = x_arr(:,2);
% end
% npts = 100;
% xi_arr = linspace(0,1,npts);
% eta_arr = xi_arr;
% x_arr = zeros(npts);
% y_arr = zeros(npts);
% f_approx = zeros(size(xi_arr));
% f_exact = zeros(size(xi_arr));
% F = zeros(size(xi_arr));
% for i = 1:length(xi_arr)
%     for j = 1:length(eta_arr)
%         xi = xi_arr(i);
%         eta = eta_arr(j);
%         [u, x, dudx, dudy] = numericalSolEval2D(xi, eta, varCol, U);
%         f_approx(i,j) = dudx - dudy + 1i*omega*(dgdx(x(1),x(2)) - dgdy(x(1),x(2)))*u;
%         f_exact(i,j) = f(x(1),x(2));
%         x_arr(i,j) = x(1);
%         y_arr(i,j) = x(2);
%         F(i,j) = u;
%     end
% end
% surf(x_arr,y_arr, 100*abs(f_approx-f_exact)/max(max(abs(f_exact))))
% % plot(xi_arr,real(F),xi_arr,real(1+exp(-1i*omega*g(xi_arr))))
% ylabel('Relative error (in %)')
% xlabel('x')
% if ~useGalerkin
%     hold on
%     plot3(cp_x,cp_y,zeros(size(cp)),'*')
% end
% d = [1, 0.8];
% d = d/norm(d);
% C = 1i*omega*(d(1)-d(2));
% c_1 = (C^2+2)/(C*(C^2+4));
% c_2 = -1/(C^2+4);
% c_3 = -1/(C^2+4);
% c_4 = 2/(C*(C^2+4));
% F_exact = @(x,y) c_1*sin(x).*cos(y) + c_2*sin(x).*sin(y) + c_3*cos(x).*cos(y) + c_4*cos(x).*sin(y);
% I_Exact2 =   integral(@(y)  F_exact(1,y).*exp(1i*omega*g(1,y)),  -1,1) ...
%            - integral(@(x) F_exact(-x,1).*exp(1i*omega*g(-x,1)), -1,1) ...
%            - integral(@(y)F_exact(-1,-y).*exp(1i*omega*g(-1,-y)),-1,1) ...
%            + integral(@(x) F_exact(x,-1).*exp(1i*omega*g(x,-1)), -1,1);
% 
% I_Exact = -2*1i*(d(1)*omega*sin(1)*cos(d(1)*omega) - cos(1)*sin(d(1)*omega))/(d(1)^2*omega^2-1)...
%            *2*(d(2)*omega*cos(1)*sin(d(2)*omega) - sin(1)*cos(d(2)*omega))/(d(2)^2*omega^2-1);
% keyboard
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I =   LevinMethod(-1,1,@(y)numericalSolEval2D(1, (y+1)/2, varCol, U), @(y)g(1,y),  @(y)dgdy(1,y),omega) ...
    - LevinMethod(-1,1,@(x)numericalSolEval2D((-x+1)/2, 1, varCol, U), @(x)g(-x,1),  @(x)-dgdx(-x,1),omega) ...
    - LevinMethod(-1,1,@(y)numericalSolEval2D(0, (-y+1)/2, varCol, U), @(y)g(-1,-y), @(y)-dgdy(-1,-y),omega) ...
    + LevinMethod(-1,1,@(x)numericalSolEval2D((x+1)/2, 0, varCol, U), @(x)g(x,-1), @(x)dgdx(x,-1),omega);
% 
% 
% function I = tempIntegrand(Y, omega, f, g)
% 
% I = zeros(size(Y));
% for i = 1:length(Y)
%     y = Y(i);
%     I(i) = f(y)*exp(1i*omega*g(y));
% end
% 






