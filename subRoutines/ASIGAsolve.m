function task = ASIGAsolve(task,omega,i_o)

switch task.misc.method
    case {'IE','IENSG','BEM','BA','ABC','MFS','PML'}        
        %% SOLVE SYSTEM
        useA = isfield(task,'A');
        t_start_sol = tic;
        if task.misc.printLog
            fprintf(['\n%-' num2str(task.misc.stringShift) 's'], ['Creating preconditioner (' task.sol.preconditioner ') ... '])
        end
        if useA
            task = createPreconditioner(task);
        end
        task.timePreconditioner = toc;
        if task.misc.printLog
            fprintf('using %12f seconds.', task.timePreconditioner)
        end

        if task.misc.printLog && ~task.rom.useROM
            fprintf(['\n%-' num2str(task.misc.stringShift) 's'], 'Solving system of equations ... ')
        end
        if useA
            task.sol.maxit = min(task.sol.maxit, size(task.A,1));
        end
        switch task.sol.solver
            case 'LU'
                if ~(strcmp(task.sol.preconditioner,'diag') || strcmp(task.sol.preconditioner,'none'))
                    error('not implemented')
                end
                if task.rom.useROM
                    task = computeROMderivatives(task);
                else
                    if useA
                        if i_o == 1 && strcmp(task.misc.formulation,'Sweep')
                            task.UU = zeros(size(task.A,1),numel(omega));
                        end
                        if strcmp(task.misc.scatteringCase,'Sweep')
                            task.UU(:,i_o) = task.Pinv*((task.Pinv*task.A*task.Pinv)\(task.Pinv*task.FF));
                        else
                            task.UU = task.Pinv*((task.Pinv*task.A*task.Pinv)\(task.Pinv*task.FF));
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
            case 'gmres'
                task.UU = zeros(size(task.FF));
                if task.sol.restart > size(task.A,1) % RESTART should be bounded by SIZE(task.A,1).
                    task.sol.restart = size(task.A,1);
                end

                for i = 1:size(task.FF,2)
                    [task.UU(:,i),task.sol.flag,task.sol.relres,task.sol.iter] = gmres(task.A,task.FF(:,i),task.sol.restart,task.sol.tol,task.sol.maxit,task.L_A,task.U_A);
%                     gmres_parf_obj = parfeval(backgroundPool,@gmres,1,task.A,task.FF(:,i),task.sol.restart,task.sol.tol,task.sol.maxit,task.L_A,task.U_A);
%                     task.UU(:,i) = fetchOutputs(gmres_parf_obj);
                end
                if task.sol.flag
                    warning('gmres did not converge to the required tolerance')
                end
            otherwise
                task.UU = zeros(size(task.FF));
                for i = 1:size(task.FF,2)
                    eval(['[task.UU(:,i),task.sol.flag,task.sol.relres,task.sol.iter] = ' task.sol.solver '(task.A,task.FF(:,i),task.sol.tol,task.sol.maxit,task.L_A,task.U_A);'])
                end
        end
        if useA
            dofs = size(task.A,1);
        else
            dofs = size(task.varCol{1}.A_K,1);
        end
        task.dofs = dofs;
        if task.misc.printLog && ~task.rom.useROM
            fprintf('using %12f seconds.', toc)
        end
        task.timeSolveSystem = toc(t_start_sol);
        if task.misc.computeCondNumber && (size(task.A,1) == size(task.A,2))
            if task.misc.printLog
                fprintf(['\n%-' num2str(task.misc.stringShift) 's'], 'Calculating condition number ... ')
            end
            rng('default') % for reproducibility in condest
            condNumber = condest(task.A);
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
            task.varCol = rmfields(task.varCol,{'A_K','A_M','A_C','Ainf','Ainf1','Ainf2'});
            task = rmfields(task,{'A','FF','Pinv','P_rightinv','P_right'});
            if ~task.rom.useROM
                task = rmfields(task,{'A0','A1','A2','A4'});
            end
        end
        % fprintf('\nMemory ratio = %f', ((fluid.degree(1)+1)^6*task.varCol{1}.noElems)/nnz(A_fluid_o))
end