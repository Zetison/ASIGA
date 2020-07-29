function I = LevinMethod(a,b,f,g,dgdx,omega)
useGalerkin = 0;
spline = getRodData(a,b-a);
degreeElev = 10;
newXiKnots = 2;

spline = elevateNURBSdegree(spline,degreeElev);
spline = insertKnotsInNURBS(spline,insertUniform2(spline.knots, newXiKnots));

varCol.dimension = 1;    
varCol = convertNURBS(spline, varCol);
varCol = generateIGA1DMesh(varCol);

p_xi = spline.degree;
noElems = varCol.noElems;
controlPts = varCol.controlPts;
noDofs = spline.number;

element = varCol.element;
weights = varCol.weights;
Xi = spline.knots;

if useGalerkin
    elRangeXi = varCol.elRangeXi;
    index = varCol.index;
    [W1D,Q1D] = gaussianQuadNURBS(p_xi+1);
    n_en = p_xi+1;
    Avalues = zeros(n_en^2,noElems); 
    spIdxRow = zeros(n_en^2,noElems);
    spIdxCol = zeros(n_en^2,noElems);
    F_indices = zeros(n_en,noElems); 
    Fvalues = zeros(n_en,noElems); 
    for e = 1:noElems
        idXi = index(e,1);

        Xi_e = elRangeXi(idXi,:);

        sctr = element(e,:);

        pts = controlPts(sctr,:);
        a_e = zeros(n_en);
        f_e = zeros(n_en,1);
        
        J_2 = 0.5*(Xi_e(2)-Xi_e(1));
        for gp = 1:size(W1D,1)
            pt = Q1D(gp,:);
            wt = W1D(gp);

            xi  = parent2ParametricSpace(Xi_e, pt(1));

            [R, dRdxi] = NURBS1DBasis2(xi, p_xi, Xi, weights);
            J_1 = dRdxi*pts;
        	dRdX = dRdxi/J_1;
            x = R*pts;
            a_e = a_e + R'*(dRdX + 1i*omega*dgdx(x)*R)* abs(J_1) * J_2 * wt;  
            f_e = f_e + R'*f(x) * abs(J_1) * J_2 * wt;
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
else
    uniqueXi = unique(Xi);

    c = aveknt(Xi, p_xi+1);

    A = zeros(noDofs);
    FF = zeros(noDofs,1);
    x_arr = zeros(noDofs,1);
    for i = 1:length(c)
        xi = c(i);
        xi_idx = findKnotSpan(noElems, 0, xi, uniqueXi);
        e = xi_idx;
        sctr = element(e,:);
        [R, dRdxi] = NURBS1DBasis2(xi, p_xi, Xi, weights);
        pts = controlPts(sctr);
        J = dRdxi*pts;
        dRdX = dRdxi/J;
        x = R*pts;
        x_arr(i) = x;
        A(i,sctr) = dRdX + 1i*omega*dgdx(x)*R;
        FF(i) = f(x);
    end

    U = A\FF;
end
v_b = numericalSolEval1D(1, varCol, U);
v_a = numericalSolEval1D(0, varCol, U);
I = v_b*exp(1i*omega*g(b)) - v_a*exp(1i*omega*g(a));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% xi_arr = linspace(0,1,1000);
% x_arr = zeros(noDofs,1);
% f_approx = zeros(size(xi_arr));
% F = zeros(size(xi_arr));
% for i = 1:length(xi_arr)
%     xi = xi_arr(i);
%     [u, x, dudx] = numericalSolEval1D(xi, varCol, U);
%     f_approx(i) = dudx + 1i*omega*dgdx(x)*u;
%     x_arr(i) = x;
%     F(i) = u;
%     
% end
% semilogy(x_arr,100*abs(f_approx.'-f(x_arr))/max(abs(f(x_arr))))
% % plot(xi_arr,real(F),xi_arr,real(1+exp(-1i*omega*g(xi_arr))))
% ylabel('Relative error (in %)')
% xlabel('x')
% keyboard
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




