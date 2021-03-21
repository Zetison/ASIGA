function p_h = calculateScatteredPressureMFS(task, P_far)

k = task.varCol{1}.k;
U = task.varCol{1}.U;

switch task.varCol{1}.formulation
    case 'SS'
        if task.ffp.plotFarField
            y = task.varCol{1}.x_0;
            R = norm2(P_far-y);
            theta = acos((P_far(:,3)-y(3))./R);
            phi = atan2(P_far(:,2)-y(2),P_far(:,1)-y(1));
            eta = cos(theta);
            N = task.ie.N;
            exp_phi = exp(1i*phi*(-N:N));
            p_h = zeros(size(P_far,1),1);
            for n = 0:N
                P = legendre(n,eta);
                for m = -n:n
                    P1 = P(abs(m)+1,:).';
                    p_h = p_h + 1i^(-n-1)/k*P1.*exp_phi(:,m+N+1);
                end
            end
        end
    case 'PS'
        y_s = task.mfs.y_s;
        if plotFarField
            x_hat = elementProd(1./norm2(P_far),P_far);
            switch task.misc.scatteringCase
                case {'BI','Sweep'}
                    p_h = 1/(4*pi)*exp(-1i*k*x_hat*y_s.')*U;
                case 'MS'
                    p_h = sum(1/(4*pi)*exp(-1i*k*x_hat*y_s.').*U.',2);
            end
        else
            n_cp = numel(U);
            y_s = task.varCol{1}.y_s;
            rs = zeros(size(P_far,1),n_cp);
            parfor j = 1:n_cp
                xmys = elementAddition(-y_s(j,:), P_far);
                rs(:,j) = norm2(xmys);
            end
            switch task.misc.scatteringCase
                case {'BI','Sweep'}
                    p_h = Phi_k(rs,k)*U;
                case 'MS'
                    p_h = sum(Phi_k(rs,k)*U.',2);
            end
        end

end
if numel(k) > 1
    p_h = p_h.';
end

