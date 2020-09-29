if varCol{1}.boundaryMethod
    if task.calculateSurfaceError
        tic 
        if ~runTasksInParallel
            fprintf(['\n%-' num2str(stringShift) 's'], 'Calculating surface error ... ')
        end
        surfaceError = calcSurfErrorBndryMethodVec(varCol, Uc, task.LpOrder);
        if ~runTasksInParallel
            fprintf('using %12f seconds.', toc)   
            fprintf('\nSurface error = %.16g', surfaceError)
        end
        if i_f == 1
            task.results.surfaceError = zeros(1,size(k,2));
        end
        task.results.surfaceError(i_f) = surfaceError;
    end
    if task.calculateSurfEnrgErr
        tic 
        if ~runTasksInParallel
            fprintf(['\n%-' num2str(stringShift) 's'], 'Calculating surface energy error ... ')
        end
        energyError = calcEnergyErrorBEM(varCol{1}, Uc{1});
        if ~runTasksInParallel
            fprintf('using %12f seconds.', toc)   
            fprintf('\nSurface energy error = %.16g', energyError)
        end
        if i_f == 1
            task.results.energyError = zeros(1,size(k,2));
        end
        task.results.energyError(i_f) = energyError;
    end
else
    if task.calculateSurfaceError
        tic 
        if ~runTasksInParallel
            fprintf(['\n%-' num2str(stringShift) 's'], 'Calculating surface error ... ')
        end
        if strcmp(varCol{1}.coreMethod,'SEM')
            surfaceError = calcSurfErrorSEM(varCol{1}, Uc{1}, task.LpOrder);
        else
            surfaceError = calcSurfErrorVec(varCol{1}, Uc{1}, task.LpOrder);
        end
        if ~runTasksInParallel
            fprintf('using %12f seconds.', toc) 
            fprintf('\nSurface error = %.16g', surfaceError)
        end
        if i_f == 1
            task.results.surfaceError = zeros(1,size(k,2));
        end
        task.results.surfaceError(i_f) = surfaceError;
    end
    if task.calculateVolumeError
        tic 
        if ~runTasksInParallel
            fprintf(['\n%-' num2str(stringShift) 's'], 'Calculating error ... ')
        end
        if strcmp(varCol{1}.coreMethod,'SEM')
            [L2Error, H1Error, H1sError, energyError] = calcErrorSEM(varCol{1}, Uc{1}(1:varCol{1}.noDofs));
        else
            [L2Error, H1Error, H1sError, energyError] = calcErrorVec(varCol, Uc);
        end
        if ~runTasksInParallel
            fprintf('using %12f seconds.', toc)
            fprintf('\nVolume L2-error = %.16g', L2Error)
            fprintf('\nVolume H1-error = %.16g', H1Error)
            fprintf('\nVolume Energy-error = %.16g', energyError)
        end
        if i_f == 1
            task.results.energyError = zeros(1,size(k,2));
            task.results.L2Error = zeros(1,size(k,2));
            task.results.H1Error = zeros(1,size(k,2));
            task.results.H1sError = zeros(1,size(k,2));
        end
        task.results.energyError(i_f) = energyError;
        task.results.L2Error(i_f) = L2Error;
        task.results.H1Error(i_f) = H1Error;
        task.results.H1sError(i_f) = H1sError;
    end  
end