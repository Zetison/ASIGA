function task = createPreconditioner(task)

task.Pinv = spdiags(1./sqrt(diag(task.A)),0,size(task.A,1),size(task.A,2));
switch task.sol.preconditioner
    case 'ilu'
        [task.L_A,task.U_A] = ilu(task.A,struct('type',task.sol.ilutype,'droptol',task.sol.droptol));
    case 'CSLP'
        beta = task.sol.beta_CSLP;  
        A_M = sparse(task.varCol{1}.noCols_tot,task.varCol{1}.noRows_tot);
        omega = task.misc.omega;
        Aindices = task.varCol{1}.Aindices;
        for i = 1:numel(task.varCol)
            if task.misc.symmetric && strcmp(task.varCol{i}.media,'solid')
                massMatrixScale = omega^4*task.varCol{i}.massMatrixScale;   
            else
                massMatrixScale = omega^2*task.varCol{i}.massMatrixScale;   
            end
            A_M(Aindices{i,1},Aindices{i,2}) = 1i*beta*abs(massMatrixScale)*task.varCol{i}.A_M;
        end
        A_M(task.varCol{1}.allDofsToRemove,:) = [];
        A_M(:,task.varCol{1}.allDofsToRemove) = [];
        [task.L_A,task.U_A] = ilu(task.A + A_M,struct('type',task.sol.ilutype,'droptol',task.sol.droptol));
%         task.Pinv = spdiags(1./sqrt(diag(task.A)),0,size(task.A,1),size(task.A,2));
%         [task.L_A,task.U_A] = ilu(task.Pinv*(task.A + A_M)*task.Pinv,struct('type',task.sol.ilutype,'droptol',task.sol.droptol));
% %         [task.L_A,task.U_A] = ilu(task.Pinv*task.A*task.Pinv,struct('type',task.sol.ilutype,'droptol',task.sol.droptol));
    case 'SSOR'
        D_SSOR = spdiags(spdiags(task.A,0),0,size(task.A,1),size(task.A,2));
        D_SSORinv = spdiags(1./spdiags(task.A,0),0,size(task.A,1),size(task.A,2));
        F_SSOR = -triu(task.A);
        E_SSOR = -tril(task.A);
        omega_SSOR = 1.5;
        task.L_A = (D_SSOR-omega_SSOR*E_SSOR)*D_SSORinv;
        task.U_A = D_SSOR-omega_SSOR*F_SSOR;
    case 'diag'
        task.Pinv = spdiags(1./sqrt(diag(task.A)),0,size(task.A,1),size(task.A,2));
    case 'none'
        task.Pinv = speye(size(task.A,1),size(task.A,2));
end