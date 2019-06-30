function tau = computeTau(varCol)
patches = varCol.patches;
dofsToRemove = varCol.dofsToRemove;
lambda = 2*pi/varCol.k;
noDofs = varCol.noDofs;
n_cp = noDofs - length(dofsToRemove);
cp = zeros(n_cp,3);
counter = 1;
counter2 = 1;
for patch = 1:varCol.noPatches
    nurbs = patches{patch}.nurbs;
    n_xi = nurbs.number(1);
    n_eta = nurbs.number(2);
    p_xi = nurbs.degree(1);
    p_eta = nurbs.degree(2);
    Xi_y = nurbs.knots{1};
    Eta_y = nurbs.knots{2};

    for j = 1:n_eta
        eta_bar = sum(Eta_y(j+1:j+p_eta))/p_eta;
        for i = 1:n_xi
            if ~any(dofsToRemove == counter)
                xi_bar = sum(Xi_y(i+1:i+p_xi))/p_xi;
                cp(counter2,:) = evaluateNURBS(nurbs,[xi_bar, eta_bar]);
                counter2 = counter2 + 1;
            end
            counter = counter + 1;
        end
    end
end

d_min = zeros(n_cp,1);
parfor i = 1:n_cp
    d_min(i) = min(norm2(repmat(cp(i,:),n_cp-1,1)-cp(setdiff(1:end,i),:)));
end
tau = lambda(1)/max(d_min);