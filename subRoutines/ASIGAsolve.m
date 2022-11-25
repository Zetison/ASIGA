function task = ASIGAsolve(task,omega,i_o)

omega_i = task.misc.omega;
switch task.misc.method
    case {'IE','IENSG','BEM','BA','ABC','MFS','PML'}
        tic
        if task.misc.printLog
            fprintf(['\n%-' num2str(task.misc.stringShift) 's'], 'Collecting matrices ... ')
        end
        useA = ~(task.noDomains == 1 && strcmp(task.misc.method,'BEM'));
        if useA
            [task,FF,A0,A1,A2,A4] = collectMatrices(task);
            if strcmp(task.misc.method,'BA')
                A = A2;
            else
                A = A0 + omega_i*A1 + omega_i^2*A2 + omega_i^4*A4;
            end
        else
            allDofsToRemove = task.varCol{1}.dofsToRemove;
            if ~(strcmp(task.misc.method,'BEM') && strcmp(task.misc.formulation(1),'C'))
                task.varCol{1}.A_K(allDofsToRemove,:) = [];
                task.varCol{1}.FF(allDofsToRemove,:) = [];
            end
            task.varCol{1}.A_K(:,allDofsToRemove) = [];
        end
        task.timeCollectingMatrices = toc;
        if task.misc.printLog
            fprintf('using %12f seconds.', task.timeCollectingMatrices)
        end
        
        if task.misc.printLog
            fprintf(['\n%-' num2str(task.misc.stringShift) 's'], 'Total time building system ... ')
            fprintf(' was  %12f seconds.', task.timeBuildSystem+task.timeCollectingMatrices)
        end
        
        %% SOLVE SYSTEM
        tic
        if task.misc.printLog
            fprintf(['\n%-' num2str(task.misc.stringShift) 's'], ['Creating preconditioner (' task.sol.preconditioner ') ... '])
        end
        if useA
            switch task.sol.preconditioner
                case 'ilu'
                    [L_A,U_A] = ilu(A,struct('type','nofill'));
                case 'SSOR'
                    D_SSOR = spdiags(spdiags(A,0),0,size(A,1),size(A,2));
                    D_SSORinv = spdiags(1./spdiags(A,0),0,size(A,1),size(A,2));
                    F_SSOR = -triu(A);
                    E_SSOR = -tril(A);
                    omega_SSOR = 1.5;
                    L_A = (D_SSOR-omega_SSOR*E_SSOR)*D_SSORinv;
                    U_A = D_SSOR-omega_SSOR*F_SSOR;
                case 'diag'
                    Pinv = spdiags(1./sqrt(diag(A)),0,size(A,2),size(A,2));
            end
        end
        if task.misc.printLog
            fprintf('using %12f seconds.', toc)
        end
        tic
        if task.misc.printLog
            fprintf(['\n%-' num2str(task.misc.stringShift) 's'], 'Solving system of equations ... ')
        end
        switch task.sol.solver
            case 'gmres'
                noRestarts = []; % change this to save memory
                % noRestarts = 2;
                [task.UU,~,~,it1,rv1] = gmres(A,FF,noRestarts,1e-20,1000,L_A,U_A);
            case 'LU'
                if ~strcmp(task.sol.preconditioner,'diag')
                    error('not implemented')
                end
                if task.rom.useROM
                    dAdomega = cell(1,4);
                    dAdomega{1} = A1 + 2*omega_i*A2 + 4*omega_i^3*A4;
                    dAdomega{2} = 2*A2 + 12*omega_i^2*A4;
                    dAdomega{3} = 24*omega_i*A4;
                    dAdomega{4} = 24*A4;
                    if task.rom.useDGP && i_o == numel(omega)
                        task.A0 = A0;
                        task.A1 = A1;
                        task.A2 = A2;
                        task.A4 = A4;
                    end
                    task.UU = zeros(size(FF));
                    dA = decomposition(Pinv*A*Pinv,'lu');
                    fprintf('using %12f seconds.', toc)
                    fprintf(['\n%-' num2str(task.misc.stringShift) 's'], 'Computing ROM solution ... ')
                    for i = 1:task.noRHSs
                        j = i-1;
                        b = FF(:,i);
                        for k = 1:min(j,numel(dAdomega))
                            b = b - nchoosek(j,k)*dAdomega{k}*task.UU(:,i-k);
                        end
                        task.UU(:,i) = Pinv*(dA\(Pinv*b));
                    end
                else
                    if useA
                        if i_o == 1 && strcmp(task.misc.formulation,'Sweep')
                            task.UU = zeros(size(A,1),numel(omega));
                        end
                        if strcmp(task.misc.scatteringCase,'Sweep')
                            task.UU(:,i_o) = Pinv*((Pinv*A*Pinv)\(Pinv*FF));
                        else
                            task.UU = Pinv*((Pinv*A*Pinv)\(Pinv*FF));
                        end
                    else
                        if i_o == 1 && strcmp(task.misc.formulation,'Sweep')
                            task.UU = zeros(size(task.varCol{1}.A_K,1),numel(omega));
                        end
                        if strcmp(task.misc.scatteringCase,'Sweep')
                            task.UU(:,i_o) = task.varCol{1}.A_K\task.varCol{1}.FF;
                        else
                            task.UU = task.varCol{1}.A_K\task.varCol{1}.FF;
                        end
                    end
                end
            otherwise
                eval(['[UU,~,~,it1,rv1] = cgs' task.sol.solver '(A,FF,1e-20,1000,L_A,U_A);'])
        end
        if useA
            dofs = size(A,1);
        else
            dofs = size(task.varCol{1}.A_K,1);
        end
        task.dofs = dofs;
        if task.misc.printLog
            fprintf('using %12f seconds.', toc)
        end
        task.timeSolveSystem = toc;
        if task.misc.computeCondNumber && (size(A,1) == size(A,2))
            if task.misc.printLog
                fprintf(['\n%-' num2str(task.misc.stringShift) 's'], 'Calculating condition number ... ')
            end
            rng('default') % for reproducibility in condest
            condNumber = condest(A);
%             condNumber = cond(full(A))
            if task.misc.printLog
                fprintf('using %12f seconds.', toc)
                fprintf('\nCondition number = %d', condNumber)
            end
        else
            condNumber = NaN;
        end
        task.results.cond_number = condNumber;

        if task.misc.printLog
            fprintf('\nNumber of degrees of freedom = %d', dofs)
        end
        if task.misc.clearGlobalMatrices
            task.varCol = rmfields(task.varCol,{'A_K','A_M','A_C','FF','Ainf','Ainf1','Ainf2'});
            clear A FF A0 A1 A2 A4 Pinv dA dAdomega b
        end
        % fprintf('\nMemory ratio = %f', ((fluid.degree(1)+1)^6*task.varCol{1}.noElems)/nnz(A_fluid_o))
end
        
if task.rom.useROM && strcmp(task.misc.scatteringCase,'Sweep')
    task.U_sweep{i_o} = task.UU;
    task = postProcessSolution(task);
end