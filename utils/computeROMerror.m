function relError = computeROMerror(task, omega_cj)
n_freq = numel(omega_cj);
relError = zeros(1,n_freq);
for i = 1:n_freq
    omega = omega_cj(i);
    task.misc.omega = omega;

    task = getAnalyticSolutions(task);
    task = buildRHS(task,true);
    task = collectMatrices(task,false,true,true);
    FFm = task.V'*task.FF;
    Am = task.A0_am + omega*task.A1_am + omega^2*task.A2_am + omega^4*task.A4_am;
    if strcmp(task.sol.preconditioner,'none')
        UU = task.V*(Am\FFm);
    else
        Pinv = spdiags(1./sqrt(diag(Am)),0,size(Am,1),size(Am,2));
        UU = task.V*(Pinv*((Pinv*Am*Pinv)\(Pinv*FFm)));
    end

    task = collectMatrices(task,false,false);
    task = createPreconditioner(task);
    UU_ref = task.Pinv*((task.Pinv*task.A*task.Pinv)\(task.Pinv*task.FF));

    relError(i) = norm(UU - UU_ref)/norm(UU_ref);
end