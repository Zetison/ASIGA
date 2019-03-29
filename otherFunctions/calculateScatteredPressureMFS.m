function p_h = calculateScatteredPressureMFS(varCol, U, P_far, plotFarField)

k = varCol.k;
y_s = varCol.y_s;

Phi_k = @(r) exp(1i*k*r)./(4*pi*r);



if plotFarField
    x_hat = elementProd(1./norm2(P_far),P_far);
    switch varCol.scatteringCase
        case {'BI','Sweep'}
            p_h = 1/(4*pi)*exp(-1i*k*x_hat*y_s.')*U;
        case 'MS'
            p_h = sum(1/(4*pi)*exp(-1i*k*x_hat*y_s.').*U.',2);
    end
else
    n_cp = numel(U);
    y_s = varCol.y_s;
    rs = zeros(size(P_far,1),n_cp);
    parfor j = 1:n_cp
        xmys = elementAddition(-y_s(j,:), P_far);
        rs(:,j) = norm2(xmys);
    end
    switch varCol.scatteringCase
        case {'BI','Sweep'}
            p_h = Phi_k(rs)*U;
        case 'MS'
            p_h = sum(Phi_k(rs)*U.',2);
    end
end


