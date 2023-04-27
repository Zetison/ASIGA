function task = Hetmaniuk2012raa_P1(task)


oldprintLog = task.misc.printLog;
if oldprintLog
    tic
    fprintf(['\n%-' num2str(task.misc.stringShift) 's'], 'Computing ROM basis adaptively ')
end
omega = task.misc.omega;
n_c = task.rom.n_c;
task.misc.printLog = false;

% Algorithm P1 in Hetmaniuk2013aas available at https://www.doi.org/10.1002/nme.4436
task.V = [];
task.U = [];
task.rom.omega_U = [];
task.rom.J_U = [];
omega_T = zeros(1,0); % Initiate with an empty row vector
omega_P = sort(task.rom.omega);
omega_T_new = [omega_P(1),omega_P(end)];
task.rom.history = struct();
J_P = [];
counter = 1;
while ~isempty(omega_T_new)
    J_new = zeros(1,numel(omega_T_new));
    for p = 1:numel(omega_T_new)
        % Algorithm P2 in Hetmaniuk2013aas
        [task, J_new(p)] = Hetmaniuk2012raa_P2(task, union(omega_T_new,omega_T), omega_T_new(p));
    end
    fprintf('\nAdded %d new snapshots to ROM basis', numel(omega_T_new))
    [omega_T,sortIdx] = sort([omega_T,omega_T_new]);
    omega_T_new = zeros(1,0); % Initiate with an empty row vector
    J_P = [J_P,J_new];
    J_P = J_P(sortIdx);
    task.rom.history(counter).residual = [];
    task.rom.history(counter).omega = [];
    for p = 1:numel(omega_T)-1
        omega_p = omega_T(p);
        omega_pp1 = omega_T(p+1);
        j = 1:n_c;
        omega_cj = omega_p + j/(n_c+1)*(omega_pp1 - omega_p);
        residual = computeROMresidual(task, omega_cj);
        task.rom.history(counter).residual = [task.rom.history(counter).residual, residual];
        task.rom.history(counter).omega = [task.rom.history(counter).omega, omega_cj];
        [~,j_max] = max(residual);
        if residual(j_max) > task.rom.tolerance
            omega_T_new = union(omega_T_new,omega_cj(j_max));
        end
    end
    task.rom.history(counter).omega_P = omega_T;
    task.rom.history(counter).J_P = J_P;
    task.rom.history(counter).omega_T_new = omega_T_new;
    counter = counter + 1;
end
if oldprintLog
    fprintf('\nTotal time for adaptive ROM basis computations %12f seconds.', toc)
end
if task.rom.computeROMresidualFine
    if oldprintLog
        fprintf(['\n%-' num2str(task.misc.stringShift) 's'], 'Computing final ROM residual/error ... ')
    end    
    omega = union(omega,task.rom.history(end).omega); % Add all values of omega used in adaptive algorithm to obtain interpolatory plots
    [residual, relError] = computeROMresidual(task, omega);
    task.rom.history(end).relError = relError;
    task.rom.history(end).residualFine = residual;
    task.rom.history(end).omegaFine = omega;
    if oldprintLog
        fprintf('using %12f seconds.', toc)
    end
end
task.V = task.P_rightinv*task.V; % Scale back to match scaling for task.U
task.misc.printLog = oldprintLog;