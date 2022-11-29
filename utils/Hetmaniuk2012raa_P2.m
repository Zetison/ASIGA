function task = Hetmaniuk2012raa_P2(task, omega, p)

J_min = task.rom.J_min;
J_max = task.rom.J_max;
deltaJ = task.rom.deltaJ;
noDofs = size(task.V_sweep{1},1);
noVecs = task.rom.noVecs;
tau = task.rom.tolerance;

k = 1:n_c;
if p == 1
    omega_k_star = omega(p) + k/(n_c+1)*(omega(p+1)-omega(p));
elseif p == numel(omega)
    omega_k_star = omega(p-1) + k/(n_c+1)*(omega(p)-omega(p-1));
else
    omega_k_star = [omega(p-1) + k/(n_c+1)*(omega(p)-omega(p-1)), ...
                    omega(p)   + k/(n_c+1)*(omega(p+1)-omega(p))];
end
V = zeros(noDofs,deltaJ);

residual_old = task.rom.history(end).residual;
task = computeROMresidual(task, omega_k_star);
residual = task.rom.history(end).residual;
c = numel(find(and(residual > tau, ...
                   residual./residual_old < 0.9^deltaJ)));

if all(residual < tau) || J + deltaJ > J_max || c == 0
    return
end
