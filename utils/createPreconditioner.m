function task = createPreconditioner(task)

tic
if task.misc.printLog
    fprintf(['\n%-' num2str(task.misc.stringShift) 's'], ['Creating preconditioner (' task.sol.preconditioner ') ... '])
end
switch task.sol.preconditioner
    case 'ilu'
        [task.L_A,task.U_A] = ilu(task.A,struct('type','nofill'));
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
if task.misc.printLog
    fprintf('using %12f seconds.', toc)
end