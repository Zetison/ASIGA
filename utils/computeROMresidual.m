function residual = computeROMresidual(task, omega_cj)
n_freq = numel(omega_cj);
residual = zeros(1,n_freq);
for i = 1:n_freq
    omega = omega_cj(i);
    task.misc.omega = omega;

    task = getAnalyticSolutions(task);
    task = buildRHS(task,true);
    task = collectMatrices(task,false,true);
    FFm = task.V'*task.FF;

    Am = task.A0_am + omega*task.A1_am + omega^2*task.A2_am + omega^4*task.A4_am;
    Pinv = spdiags(1./sqrt(diag(Am)),0,size(Am,1),size(Am,2));
    UU = task.V*(Pinv*((Pinv*Am*Pinv)\(Pinv*FFm)));

    LHS = (task.A0 + omega*task.A1 + omega^2*task.A2 + omega^4*task.A4)*UU;
    residual(i) = norm(LHS - task.FF)/norm(task.FF);
end