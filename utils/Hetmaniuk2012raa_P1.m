function task = Hetmaniuk2012raa_P1(task)


oldprintLog = task.misc.printLog;
if oldprintLog
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
task.timeBuildROM = 0;
while ~isempty(omega_T_new)
    t_Hetmaniuk2012raa_P2_start = tic;
    J_new = zeros(1,numel(omega_T_new));
    for p = 1:numel(omega_T_new)
        % Algorithm P2 in Hetmaniuk2013aas
        [task, J_new(p)] = Hetmaniuk2012raa_P2(task, union(omega_T_new,omega_T), omega_T_new(p));
    end
    t_toc = toc(t_Hetmaniuk2012raa_P2_start);
    fprintf('\nAdded %d new snapshots to ROM basis using %12f seconds.', numel(omega_T_new), t_toc)
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
    task.timeBuildROM = task.timeBuildROM + t_toc;
%     plotROMresiduals(task);
end
if oldprintLog
    fprintf('\nTotal time for adaptive ROM basis computations: %12f seconds.', task.timeBuildROM)
end
if task.rom.computeROMresidualFine
    tic
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

if strcmp(task.sol.preconditioner,'CSLP')
    for i = 1:numel(task.varCol)
        task.varCol{i} = rmfield(task.varCol{i},'A_M');
    end
end

%% Print history to file
history = task.rom.history;
options.xlabel = 'k';
options.ylabel = 'residual';
c_f = task.varCol{1}.c_f;
for i = 1:numel(history)
    omega = task.rom.history(i).omega;
    residual = task.rom.history(i).residual;

    omega_P = history(i).omega_P;
    J_P = history(i).J_P;
    options.x = omega_P.'/c_f;
    options.y = 1e-15*ones(size(omega_P)).';
    printResultsToFile([task.resultsFolder, '/', task.saveName '_InterpolationPoints_iter' num2str(i)], options)

    options.x = omega_P.'/c_f;
    options.y = J_P.';
    options.ylabel = 'J_P';
    printResultsToFile([task.resultsFolder, '/', task.saveName '_noDerivatives_iter' num2str(i)], options)

    options.x = [omega_P(1),omega_P(end)].'/c_f;
    options.y = task.rom.tolerance*ones(2,1);
    options.ylabel = 'residual';
    printResultsToFile([task.resultsFolder, '/', task.saveName '_tolerance'], options)

    options.x = omega.'/c_f;
    options.y = residual.';
    printResultsToFile([task.resultsFolder, '/', task.saveName '_residual_iter' num2str(i)], options)

    omega_T_new = history(i).omega_T_new;
    [~,idx] = ismember(omega,omega_T_new);

    options.x = omega_T_new.'/c_f;
    options.y = residual(logical(idx)).';
    printResultsToFile([task.resultsFolder, '/', task.saveName '_NewInterpolationPoints_iter' num2str(i)], options)
end
if isfield(task.rom.history(end),'residualFine')
    residual = task.rom.history(end).residualFine;
    omegaFine = task.rom.history(end).omegaFine;

    options.x = omegaFine.'/c_f;
    options.y = residual.';
    printResultsToFile([task.resultsFolder, '/', task.saveName '_FinalRelativeResidual'], options)
end
if isfield(task.rom.history(end),'relError')
    relError = task.rom.history(end).relError;
    options.x = omegaFine.'/c_f;
    options.y = relError.';
    options.xlabel = 'k';
    options.ylabel = 'error';
    printResultsToFile([task.resultsFolder, '/', task.saveName '_FinalRelativeError'], options)
end