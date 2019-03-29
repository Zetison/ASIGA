
if boundaryMethod
    if calculateSurfaceError
        tic 
        if ~runTasksInParallel
            fprintf(['\n%-' num2str(stringShift) 's'], 'Calculating surface error ... ')
        end
        surfaceError = calcSurfErrorBndryMethod(varCol, U_fluid_o, task.LpOrder);
        if ~runTasksInParallel
            fprintf('using %12f seconds.', toc)   
            fprintf('\nSurface error = %g', surfaceError)
        end
        task.results.surfaceError = surfaceError;
    end
else
    if calculateSurfaceError
        tic 
        if ~runTasksInParallel
            fprintf(['\n%-' num2str(stringShift) 's'], 'Calculating surface error ... ')
        end
        surfaceError = calcSurfErrorVec(varCol, U_fluid_o, task.LpOrder);
        if ~runTasksInParallel
            fprintf('using %12f seconds.', toc) 
            fprintf('\nSurface error = %g', surfaceError)
        end
        task.results.surfaceError = surfaceError;
    end
    if calculateVolumeError
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
            [L2Error, H1Error, H1sError, energyError] = calcError(varColCell, U_cell, e3Dss_options);
        end
        if ~runTasksInParallel
            fprintf('using %12f seconds.', toc)
            fprintf('\nVolume L2-error = %f', L2Error)
            fprintf('\nVolume H1-error = %f', H1Error)
            fprintf('\nVolume Energy-error = %f', energyError)
        end
        task.results.energyError = energyError;
        task.results.L2Error = L2Error;
        task.results.H1Error = H1Error;
        task.results.H1sError = H1sError;
    end  
end