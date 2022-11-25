function task = ASIGAassembly(task,t_start)


omega_i = task.misc.omega;
task = getAnalyticSolutions(task);
switch task.misc.method
    case {'IE','ABC','IENSG','PML'}
        task.timeBuildSystem = 0;
        for i_domain = 1:task.noDomains
            if ~(strcmp(task.misc.method,'IENSG') && task.iem.boundaryMethod && i_domain == 1)
                tic          
                if task.misc.printLog
                    fprintf(['\n%-' num2str(task.misc.stringShift) 's'], ['Building matrices for domain ' num2str(i_domain) ' ... '])
                end
                if strcmp(task.misc.coreMethod,'SEM')  
                    if task.noDomains > 1
                        error('Not implemented')
                    end  
                    task = buildSEMMatrices(task);
                else  
                    task = buildMatrices(task,i_domain);
                end
                task.timeBuildSystem = task.timeBuildSystem + toc;
                if task.misc.printLog
                    fprintf('using %12f seconds.', toc)
                end
            end
        end
        if ~strcmp(task.misc.coreMethod,'SEM')
            if strcmp(task.misc.method,'IE')
                % Add contribution from infinite elements
                tic
                if task.misc.printLog
                    fprintf(['\n%-' num2str(task.misc.stringShift) 's'], 'Building infinite element matrix ... ')
                end
                task = buildIEmatrix(task);
                
                task.timeBuildSystem = task.timeBuildSystem + toc;
                if task.misc.printLog
                    fprintf('using %12f seconds.', toc)
                end
            elseif strcmp(task.misc.method,'IENSG')
                tic
                if task.misc.printLog
                    fprintf(['\n%-' num2str(task.misc.stringShift) 's'], 'Building infinite element matrix ... ')
                end
                chimax = task.varCol{1}.chimax;
                chimin = task.varCol{1}.chimin;
                if abs(chimax - chimin)/abs(chimax) < 100*eps % Exploit tensor product structure of IEM
                    task.misc.r_a = mean([chimax,chimin]);
                    task = buildIEmatrix(task);
                else
                    task = infElementsNonSepGeom(task);  
                end
                if task.misc.printLog
                    fprintf('using %12f seconds.', toc)
                end
            elseif strcmp(task.misc.method,'ABC')
                % Add contribution from infinite elements
                tic
                if task.misc.printLog
                    fprintf(['\n%-' num2str(task.misc.stringShift) 's'], 'Building ABC matrix ... ')
                end
                task.varCol{1} = addABC(task.varCol{1}); 

                task.timeBuildSystem = task.timeBuildSystem + toc;
                if task.misc.printLog
                    fprintf('using %12f seconds.', toc)
                end
            end  
            % Apply coupling conditions  
            if task.noDomains > 1 
                tic
                if task.misc.printLog
                    fprintf(['\n%-' num2str(task.misc.stringShift) 's'], 'Building coupling matrix ... ')
                end
                for i_domain = 2:task.noDomains
                    task = applyCouplingConditionPatches(task,i_domain);
                end
                if task.misc.printLog
                    fprintf('using %12f seconds.', toc)
                end
                task.timeBuildSystem = task.timeBuildSystem + toc;
            end

            % Apply Neumann conditions
            tic
            if task.misc.printLog
                fprintf(['\n%-' num2str(task.misc.stringShift) 's'], 'Building RHS vector ... ')
            end
            for i_domain = 1:min(task.noDomains,2)
                task = applyNeumannCondition(task,i_domain);
                if task.misc.symmetric && ~task.rom.useROM && strcmp(task.varCol{i_domain}.media,'solid')
                    task.varCol{i_domain}.FF = omega_i^2*task.varCol{i_domain}.FF;
                end
            end
            if task.misc.printLog
                fprintf('using %12f seconds.', toc)
            end
            task.timeBuildSystem = task.timeBuildSystem + toc;
        end
    case 'BEM'
        tic
        if task.misc.printLog
            fprintf(['\n%-' num2str(task.misc.stringShift) 's'], 'Building BEM matrix ... ')
        end
        switch task.misc.formulation(1)
            case 'G' % Galerkin
                task = buildGBEMmatrix(task,1);  
            case 'C' % Collocation
                task = buildCBEMmatrix(task);  
        end
        if task.misc.printLog
            fprintf('using %12f seconds.', toc)
        end
        if strcmp(task.misc.formulation(end),'C')
            tic
            if task.misc.printLog
                fprintf(['\n%-' num2str(task.misc.stringShift) 's'], 'Building CHIEF matrix ... ')
            end
            task = buildCHIEFmatrix(task);
            if task.misc.printLog
                fprintf('using %12f seconds.', toc)
            end
        end
        task.timeBuildSystem = toc(t_start);
        if task.noDomains > 1 
            tic
            % Solid matrices
            if task.misc.printLog
                fprintf(['\n%-' num2str(task.misc.stringShift) 's'], 'Building solid matrix ... ')
            end
            task = buildMatrices(task,2);
            
            task.timeBuildSystem = task.timeBuildSystem + toc;
            if task.misc.printLog
                fprintf('using %12f seconds.', toc)
            end
        end
        if task.noDomains > 2
            if task.misc.printLog
                fprintf(['\n%-' num2str(task.misc.stringShift) 's'], 'Building inner fluid matrix ... ')
            end
            
            switch formulation(1)
                case 'G' % Galerkin
                    task = buildGBEMmatrix(task,3);  
                case 'C' % Collocation
                    error('Not implemented')
            end
            if task.misc.printLog
                fprintf('using %12f seconds.', toc)
            end
        end

        % Apply coupling conditions  
        if task.noDomains > 1 
            tic
            if task.misc.printLog
                fprintf(['\n%-' num2str(task.misc.stringShift) 's'], 'Building coupling matrix ... ')
            end  
            for i = 2:task.noDomains
                if strcmp(task.varCol{i}.media,'fluid')
                    task.varCol{i} = applyCouplingCondition_FEBE(task.varCol{i});
                end
            end
            if task.misc.printLog
                fprintf('using %12f seconds.', toc)
            end
            task.timeBuildSystem = task.timeBuildSystem + toc;
        end  
    case 'BA'
        task.timeBuildSystem = 0;
        for i = 1:task.noDomains
            tic          
            if task.misc.printLog
                fprintf(['\n%-' num2str(task.misc.stringShift) 's'], ['Building matrices for domain ' num2str(i) ' ... '])
            end
            task = buildMatrices(task,i);
            task = buildBAmatrix(task,i);
            task.timeBuildSystem = task.timeBuildSystem + toc;
            if task.misc.printLog
                fprintf('using %12f seconds.', toc)
            end
        end
%                 if task.noDomains > 1 
%                     warning('BA is not implemented in such a way that the displacement and pressure conditions at the interfaces is satisfied')
%                 end
        task.timeBuildSystem = toc;
    case 'MFS'
        tic
        if task.misc.printLog
            fprintf(['\n%-' num2str(task.misc.stringShift) 's'], 'Building MFS matrix ... ')
        end
        task = buildMFSmatrix(task);  
        task.timeBuildSystem = toc;
        if task.misc.printLog
            fprintf('using %12f seconds.', task.timeBuildSystem)
        end
end