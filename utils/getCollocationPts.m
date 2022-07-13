function [cp_p, patchIdx, n_cp] = getCollocationPts(patches,noDofs,dofsToRemove,col_C0,colMethod)

noPatches = numel(patches);
n_cp = noDofs - length(dofsToRemove);
degree = patches{1}.nurbs.degree; % assume degree is equal in all patches
if all(degree == 1)
    eps_greville = col_C0./(2*degree);
else
    eps_greville = col_C0./degree;
end
counter2 = 1;
counter = 1;
cp_p = zeros(n_cp,2);
patchIdx = zeros(n_cp,1);

for patch = 1:noPatches
    nurbs = patches{patch}.nurbs;
    degree = nurbs.degree;
    n_xi = nurbs.number(1);
    n_eta = nurbs.number(2);
    Xi_y = nurbs.knots{1};
    Eta_y = nurbs.knots{2};
    
    switch colMethod
        case 'Grev'
            for j = 1:n_eta
                eta_bar = sum(Eta_y(j+1:j+degree(2)))/degree(2);
                if Eta_y(j+1) == Eta_y(j+degree(2))
                    if Eta_y(j+1) == Eta_y(j+degree(2)+1)
                        eta_bar = eta_bar - eps_greville(2)*(Eta_y(j+degree(2)+1)-Eta_y(j));
                    else
                        eta_bar = eta_bar + eps_greville(2)*(Eta_y(j+degree(2)+1)-Eta_y(j+1));
                    end
                end
                for i = 1:n_xi
                    if ~any(dofsToRemove == counter)
                        xi_bar = sum(Xi_y(i+1:i+degree(1)))/degree(1);
                        if Xi_y(i+1) == Xi_y(i+degree(1))
                            if Xi_y(i+1) == Xi_y(i+degree(1)+1)
                                xi_bar = xi_bar - eps_greville(1)*(Xi_y(i+degree(1)+1)-Xi_y(i));
                            else
                                xi_bar = xi_bar + eps_greville(1)*(Xi_y(i+degree(1)+1)-Xi_y(i+1));
                            end
                        end

                        cp_p(counter2,:) = [xi_bar, eta_bar];
                        patchIdx(counter2) = patch;
                        counter2 = counter2 + 1;
                    end
                    counter = counter + 1;
                end
            end
        case 'CG'
            cg_xi = CauchyGalerkin(degree(1), n_xi, Xi_y);
            cg_eta = CauchyGalerkin(degree(2), n_eta, Eta_y);
            for j = 1:n_eta
                eta_bar = cg_eta(j);
                for i = 1:n_xi
                    if ~any(dofsToRemove == counter)
                        xi_bar = cg_xi(i);
                        cp_p(counter2,:) = [xi_bar, eta_bar];
                        patchIdx(counter2) = patch;
                        counter2 = counter2 + 1;
                    end
                    counter = counter + 1;
                end
            end
        case 'GL'
            cg_xi = splinesGL(Xi_y,degree(1));
            cg_eta = splinesGL(Eta_y,degree(2));
            for j = 1:n_eta
                eta_bar = cg_eta(j);
                for i = 1:n_xi
                    if ~any(dofsToRemove == counter)
                        xi_bar = cg_xi(i);
                        cp_p(counter2,:) = [xi_bar, eta_bar];
                        patchIdx(counter2) = patch;
                        counter2 = counter2 + 1;
                    end
                    counter = counter + 1;
                end
            end
    end
end