function residual = computeROMresidual(task, omega_cj)


n_freq = length(omega_cj);
batchSize = 20;
reminder = mod(n_freq,batchSize);
n_batches = (n_freq-reminder)/batchSize;
omega_cell = mat2cell(reshape(omega_cj(1:end-reminder),batchSize,[]), batchSize,ones(1,n_batches)).';
omega_cell{end+1,1} = omega_cj(end-reminder+1:end).';
n_batches = n_batches + 1;
residual = cell(size(omega_cell));
for i = 1:n_batches
    omega = omega_cell{i}.';
    task.misc.omega = omega;

    task = getAnalyticSolutions(task);
    task = buildRHS(task,true,numel(omega));
    task = collectMatrices(task,false,true,false);
    FFm = task.V'*task.FF;
    LHS = zeros(size(task.FF));
    for j = 1:numel(omega)
        task.misc.omega = omega(j);
        task = collectMatrices(task,false,false,true);
        Am = task.A0_am + omega(j)*task.A1_am + omega(j)^2*task.A2_am + omega(j)^4*task.A4_am;
        if strcmp(task.sol.preconditioner,'none')
            UU = task.V*(Am\FFm(:,j));
        else
            Pinv = spdiags(1./sqrt(diag(Am)),0,size(Am,1),size(Am,2));
            UU = task.V*(Pinv*((Pinv*Am*Pinv)\(Pinv*FFm(:,j))));
        end
    
        LHS(:,j) = task.A*UU;
    end
    residual{i} = (vecnorm(LHS - task.FF,2,1)./vecnorm(task.FF,2,1)).';
end
residual = cell2mat(residual).';