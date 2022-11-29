function task = computeROMresidual(task, omega_cj)
n_c = numel(omega_cj);
task.rom.history(end+1).residual = zeros(n_c,1);
for i = 1:n_c
    omega_i = omega_cj(i);

    task = buildRHS(task,overrideUseROM);
    [task,FF] = collectMatrices(task);
    FFm = task.V'*FF;

    Am = task.A0_am + omega_i*task.A1_am + omega_i^2*task.A2_am + omega_i^4*task.A4_am;
    Pinv = spdiags(1./sqrt(diag(Am)),0,size(Am,1),size(Am,2));
    UU = V*(Pinv*((Pinv*Am*Pinv)\(Pinv*FFm)));

    task.rom.history(end).residual(i) = norm((task.A0 + omega_i*task.A1 + omega_i^2*task.A2 + omega_i^4*task.A4)*UU - FF)/norm(FF);
end