function tau = computeTau(task)
dofsToRemove = task.varCol{1}.dofsToRemove;
k = task.misc.omega/task.varCol{1}.c_f;
lambda = 2*pi/k;
noDofs = task.varCol{1}.noDofs;
dofsToRemove(dofsToRemove > noDofs) = [];
n_cp = noDofs - length(dofsToRemove);
cp = zeros(n_cp,3);
counter = 1;
counter2 = 1;
if task.varCol{1}.boundaryMethod
    for patch = 1:numel(task.varCol{1}.nurbs)
        nurbs = task.varCol{1}.nurbs{patch};
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
else
    for patch = 1:numel(task.varCol{1}.nurbs)
        nurbs = task.varCol{1}.nurbs{patch};
        n_xi = nurbs.number(1);
        n_eta = nurbs.number(2);
        n_zeta = nurbs.number(3);
        p_xi = nurbs.degree(1);
        p_eta = nurbs.degree(2);
        p_zeta = nurbs.degree(3);
        Xi_y = nurbs.knots{1};
        Eta_y = nurbs.knots{2};
        Zeta_y = nurbs.knots{3};

        for l = 1:n_zeta
            zeta_bar = sum(Zeta_y(l+1:l+p_zeta))/p_zeta;
            for j = 1:n_eta
                eta_bar = sum(Eta_y(j+1:j+p_eta))/p_eta;
                for i = 1:n_xi
                    if ~any(dofsToRemove == counter)
                        xi_bar = sum(Xi_y(i+1:i+p_xi))/p_xi;
                        cp(counter2,:) = evaluateNURBS(nurbs,[xi_bar, eta_bar, zeta_bar]);
                        counter2 = counter2 + 1;
                    end
                    counter = counter + 1;
                end
            end
        end
    end
end

d_min = zeros(n_cp,1);
for i = 1:n_cp
    d_min(i) = min(norm2(repmat(cp(i,:),n_cp-1,1)-cp(setdiff(1:end,i),:)));
end
tau = lambda/max(d_min);