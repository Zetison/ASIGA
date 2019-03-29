function [cg, cg0] = findCGpoints2D(varCol, U)
p_xi = varCol.nurbs.degree(1);
p_eta = varCol.nurbs.degree(2);
n_xi = varCol.nurbs.number(1);
n_eta = varCol.nurbs.number(2);
Xi = varCol.nurbs.knots{1};
Eta = varCol.nurbs.knots{2};
uniqueXi = unique(Xi);
uniqueEta = unique(Eta);
k = varCol.k;

cg_xi = CauchyGalerkin(p_xi, n_xi, Xi);
cg_eta = CauchyGalerkin(p_eta, n_eta, Eta);
cg = zeros(n_xi,n_eta,2);
cg0 = zeros(n_xi,n_eta,2);

parfor j = 1:length(cg_eta)
    
    eta = cg_eta(j);
    jj = 1;
    while eta > uniqueEta(jj+1)
        jj = jj + 1;
    end
    ii = 1;
    cg_temp = zeros(length(cg_xi), 2);
    cg0_temp = zeros(length(cg_xi), 2);
    for i = 1:length(cg_xi)
        xi = cg_xi(i);
        while xi > uniqueXi(ii+1)
            ii = ii + 1;
        end
        options = optimset('TolX',1e-14);
        f = @(c) abs(R_h(c, U, varCol, k));
        cg_temp(i, :) = fminsearchbnd(f, [xi, eta],[uniqueXi(ii), uniqueEta(jj)], [uniqueXi(ii+1), uniqueEta(jj+1)], options);
        cg0_temp(i, :) = [xi, eta];
%         h = max([uniqueXi(ii+1) - uniqueXi(ii), uniqueEta(jj+1) - uniqueEta(jj)]);
%         while length(uniquetol2(cg_temp(1:i,:),1e-1*h)) < i
%             xi = uniqueXi(ii) + rand()*(uniqueXi(ii+1) - uniqueXi(ii));
%             eta = uniqueEta(jj) + rand()*(uniqueEta(jj+1) - uniqueEta(jj)); 
%             cg_temp(i, :) = fminsearchbnd(f, [xi, eta],[uniqueXi(ii), uniqueEta(jj)], [uniqueXi(ii+1), uniqueEta(jj+1)], options);
%         end
    end
    cg(:,j,:) = cg_temp;
    cg0(:,j,:) = cg0_temp;
end

function R = R_h(c, U, varCol, k)

[u, v, ~, ~, d2udx2, d2udy2] = numericalSolEval2D(c(1), c(2), varCol, U);
R = d2udx2 + d2udy2 + k^2*u + varCol.f(v);

% 
% nppe = 2*p_xi-1;
% g0_arr = zeros(noElems,nppe);
% for e = 1:noElems
%     idXi = index(e);
% 
%     Xi_e = elRangeXi(idXi,:);
% 
%     g0_arr(e,:) = linspace2(Xi_e(1),Xi_e(2),nppe);
% end  
% g0_arr = reshape(g0_arr,noElems*nppe,1);
% 
% gk_arr = zeros(size(g0_arr));
% parfor i = 1:length(g0_arr)
%     gk = g0_arr(i);
%     gk_prev = inf;
%     counter = 1;
%     while abs(gk_prev-gk) > 10*eps && counter < 200
%         gk_prev = gk;
%         if gk <= 0
%             gk = eps;
%         elseif gk > 1
%             gk = 1-eps;
%         end
%         [u, v, dudx, d2udx2, d3udx3] = numericalSolEval1D(gk, varCol, U);
%         R_h = d2udx2 + k^2*u + varCol.f(v);
%         dR_hdx = d3udx3 + k^2*dudx + varCol.dfdx(v);
%         
%         gk = gk - R_h/dR_hdx;
%         counter = counter + 1;
%     end
%     gk_arr(i) = gk;
% end
% gk_arr = sort(uniquetol(gk_arr,1e-3/noElems/nppe));