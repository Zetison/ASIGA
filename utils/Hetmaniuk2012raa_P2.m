function [task, J] = Hetmaniuk2012raa_P2(task, omega, omega_p)

J_min = task.rom.J_min;
J_max = task.rom.J_max;
deltaJ = task.rom.deltaJ;
tau = task.rom.tolerance;
n_c = task.rom.n_c;
p = find(omega == omega_p);
k = 1:n_c;
if p == 1
    omega_k_star = omega(p) + k/(n_c+1)*(omega(p+1)-omega(p));
elseif p == numel(omega)
    omega_k_star = omega(p-1) + k/(n_c+1)*(omega(p)-omega(p-1));
else
    omega_k_star = [omega(p-1) + k/(n_c+1)*(omega(p)-omega(p-1)), ...
                    omega(p)   + k/(n_c+1)*(omega(p+1)-omega(p))];
end
task.misc.omega = omega(p);
task = collectMatrices(task,false,false,true);
task = createPreconditioner(task);

noColsInV = size(task.V,2);
J = min(max(J_min,deltaJ), J_max);
residual_old = -1;
task.noRHSs = J;
shiftROM = 0;
while true 
    % Add to V deltaJ vectors computed at omega(p) by using DGP method and
    % update the interpolatory ROM
    task.misc.omega = omega(p);
    task = getAnalyticSolutions(task,shiftROM);
    task = buildRHS(task);
    task = collectMatrices(task,false,true,false);
    task = computeROMderivatives(task,shiftROM);
    task = buildDGP(task,noColsInV+shiftROM+1);

    residual = computeROMresidual(task, omega_k_star);
    if residual_old == -1
        c = -1;
    else
        c = numel(find(and(residual_old > tau, ...
                           residual./residual_old < 0.9^deltaJ)));
    end
    
    if all(residual < tau) || J + deltaJ > J_max || c == 0
        return
    end
    shiftROM = J;
    residual_old = residual;
    J = J + deltaJ;
    task.noRHSs = deltaJ;
end
