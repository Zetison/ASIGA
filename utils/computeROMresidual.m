function [residual, relError] = computeROMresidual(task, omega_cj)

n_freq = length(omega_cj);
batchSize = 20;
reminder = mod(n_freq,batchSize);
n_batches = (n_freq-reminder)/batchSize;
omega_cell = mat2cell(reshape(omega_cj(1:end-reminder),batchSize,[]), batchSize,ones(1,n_batches)).';
omega_cell{end+1,1} = omega_cj(end-reminder+1:end).';
n_batches = n_batches + 1;
computeError = task.rom.computeROMerror && nargout == 2;
if computeError
    relError = cell(size(omega_cell));
else
    relError = NaN;
end
residual = cell(size(omega_cell));
for i = 1:n_batches
    omega = omega_cell{i}.';
    task.misc.omega = omega;

    task = getAnalyticSolutions(task);
    task = buildRHS(task,true,numel(omega));
    task = collectMatrices(task,false,true,false);
    FFm = (task.P_rightinv*task.V)'*task.FF;
    UU = zeros(size(task.FF));
    UU_ref = zeros(size(task.FF));
    LHS = zeros(size(task.FF));
    for j = 1:numel(omega)
        task.misc.omega = omega(j);
        task = collectMatrices(task,false,false,true);
        Am = task.A0_am + omega(j)*task.A1_am + omega(j)^2*task.A2_am + omega(j)^4*task.A4_am;
        if strcmp(task.sol.preconditioner,'none')
            UU(:,j) = (task.P_rightinv*task.V)*(Am\FFm(:,j));
        else
            Pinv_m = spdiags(1./sqrt(diag(Am)),0,size(Am,1),size(Am,2));
            UU(:,j) = (task.P_rightinv*task.V)*(Pinv_m*((Pinv_m*Am*Pinv_m)\(Pinv_m*FFm(:,j))));
        end
        if computeError
            task = createPreconditioner(task);
            UU_ref(:,j) = task.Pinv*((task.Pinv*task.A*task.Pinv)\(task.Pinv*task.FF(:,j)));
        end
        LHS(:,j) = task.A*UU(:,j);
    end
    if computeError
        relError{i} = (vecnorm(UU - UU_ref,2,1)./vecnorm(UU_ref,2,1)).';
    end
    residual{i} = (vecnorm(LHS - task.FF,2,1)./vecnorm(task.FF,2,1)).';
end
if computeError
    relError = cell2mat(relError).';
end
residual = cell2mat(residual).';