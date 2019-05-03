
if varCol.boundaryMethod
    if task.calculateSurfaceError
        tic 
        if ~runTasksInParallel
            fprintf(['\n%-' num2str(stringShift) 's'], 'Calculating surface error ... ')
        end
        surfaceError = calcSurfErrorBndryMethodVec(varCol, U_fluid_o, task.LpOrder);
        if ~runTasksInParallel
            fprintf('using %12f seconds.', toc)   
            fprintf('\nSurface error = %.16g', surfaceError)
        end
        if i_k == 1
            task.results.surfaceError = zeros(1,size(k,2));
        end
        task.results.surfaceError(i_k) = surfaceError;
    end
else
    if task.calculateSurfaceError
        tic 
        if ~runTasksInParallel
            fprintf(['\n%-' num2str(stringShift) 's'], 'Calculating surface error ... ')
        end
        surfaceError = calcSurfErrorVec(varCol, U_fluid_o, task.LpOrder);
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
        varColCell = {varCol};
        U_cell = {U_fluid_o};
        if useSolidDomain
            varColCell{2} = varCol_solid;
            U_cell{2} = U_solid;
        end
        if useInnerFluidDomain
            varColCell{3} = varCol_fluid_i;
            U_cell{3} = U_fluid_i;
        end
        if ~runTasksInParallel
            fprintf(['\n%-' num2str(stringShift) 's'], 'Calculating error ... ')
        end
        if strcmp(varCol.coreMethod,'SEM')
            [L2Error, H1Error, H1sError, energyError] = calcErrorSEM(varCol, U_fluid_o(1:varCol.noDofs), e3Dss_options);
%                     [L2Error, H1Error, H1sError, energyError] = calcError_testSEM({varCol}, {U_fluid_o(1:varCol.noDofs)}, e3Dss_options);
        else
            [L2Error, H1Error, H1sError, energyError] = calcErrorVec(varColCell, U_cell, e3Dss_options);
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