if varCol{1}.boundaryMethod
    if task.calculateSurfaceError
        tic 
        if ~runTasksInParallel
            fprintf(['\n%-' num2str(stringShift) 's'], 'Calculating surface error ... ')
        end
        surfaceError = calcSurfErrorBndryMethodVec(varCol{1}, Uc{1}, task.LpOrder);
        if ~runTasksInParallel
            fprintf('using %12f seconds.', toc)   
            fprintf('\nSurface error = %.16g', surfaceError)
        end
        if i_k == 1
            task.results.surfaceError = zeros(1,size(k,2));
        end
        task.results.surfaceError(i_k) = surfaceError;
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
        if i_k == 1
            task.results.energyError = zeros(1,size(k,2));
        end
        task.results.energyError(i_k) = energyError;
    end
else
    if task.calculateSurfaceError
        tic 
        if ~runTasksInParallel
            fprintf(['\n%-' num2str(stringShift) 's'], 'Calculating surface error ... ')
        end
        if strcmp(varCol{1}.coreMethod,'SEM')
            surfaceError = calcSurfErrorSEM(varCol{1}, Uc{1}, task.LpOrder);
            [L2Error, H1Error, H1sError, energyError] = calcErrorSEM(varCol{1}, Uc{1}(1:varCol{1}.noDofs), e3Dss_options);
%             [L2Error, H1Error, H1sError, energyError] = calcError_testSEM({varCol}, {Uc{1}(1:varCol.noDofs)}, e3Dss_options);
        else
            surfaceError = calcSurfErrorVec(varCol{1}, Uc{1}, task.LpOrder);
        end
        if ~runTasksInParallel
            fprintf('using %12f seconds.', toc) 
            fprintf('\nSurface error = %.16g', surfaceError)
        end
        if i_k == 1
            task.results.surfaceError = zeros(1,size(k,2));
        end
        task.results.surfaceError(i_k) = surfaceError;
    end
    if task.calculateVolumeError
        tic 
        if ~runTasksInParallel
            fprintf(['\n%-' num2str(stringShift) 's'], 'Calculating error ... ')
        end
        if strcmp(varCol{1}.coreMethod,'SEM')
            [L2Error, H1Error, H1sError, energyError] = calcErrorSEM(varCol{1}, Uc{1}(1:varCol{1}.noDofs), e3Dss_options);
%                     [L2Error, H1Error, H1sError, energyError] = calcError_testSEM({varCol}, {Uc{1}(1:varCol.noDofs)}, e3Dss_options);
        else
            [L2Error, H1Error, H1sError, energyError] = calcErrorVec(varCol, Uc, e3Dss_options);
        end
        if ~runTasksInParallel
            fprintf('using %12f seconds.', toc)
            fprintf('\nVolume L2-error = %.16g', L2Error)
            fprintf('\nVolume H1-error = %.16g', H1Error)
            fprintf('\nVolume Energy-error = %.16g', energyError)
        end
        if i_k == 1
            task.results.energyError = zeros(1,size(k,2));
            task.results.L2Error = zeros(1,size(k,2));
            task.results.H1Error = zeros(1,size(k,2));
            task.results.H1sError = zeros(1,size(k,2));
        end
        task.results.energyError(i_k) = energyError;
        task.results.L2Error(i_k) = L2Error;
        task.results.H1Error(i_k) = H1Error;
        task.results.H1sError(i_k) = H1sError;
    end  
end