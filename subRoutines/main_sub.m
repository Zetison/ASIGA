function task = main_sub(task,loopParameters,runTasksInParallel,resultsFolder)

stringShift = 40;
[varCol,task] = extractTaskData(task);
task.saveName = defineFolderAndFilenames(task,loopParameters);
task.resultsFolder = resultsFolder;
if ~runTasksInParallel
    fprintf('\n\nRunning the case: %s', task.saveName)
end

varCol = setAndDefineParameters(varCol,task);

if ~runTasksInParallel
    fprintf(['\n%-' num2str(stringShift) 's'], 'Extracting CAD data ... ')
    tic
end
varCol = createNURBSmesh(varCol, task.model, task.M, task.degree);
varCol = collectVariables(varCol,task);
noDomains = numel(varCol);
if ~runTasksInParallel
    fprintf('using %12f seconds.', toc)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot geometry
% keyboard
if task.prePlot.plot3Dgeometry || task.prePlot.plot2Dgeometry
    tic
    fprintf(['\n%-' num2str(stringShift) 's'], 'Plotting geometry ... ')
    plotMeshAndGeometry(varCol,task);
    fprintf('using %12f seconds.', toc)
end

%% Build connectivity
tic
if ~runTasksInParallel
    fprintf(['\n%-' num2str(stringShift) 's'], 'Generating IGA mesh ... ')
end
for i = 1:noDomains
    varCol{i} = generateIGAmesh(varCol{i});
end

%% Find nodes to remove due to gluing and homogeneous dirichlet conditions
totNoElems = 0;
for i = 1:noDomains
    varCol{i} = findDofsToRemove(varCol{i});
    totNoElems = totNoElems + varCol{i}.noElems;
end

if ~runTasksInParallel
    fprintf('using %12f seconds.', toc)
    fprintf('\nTotal number of elements = %d', totNoElems)
end
varCol{1}.totNoElems = totNoElems;
if (task.prePlot.plot3Dgeometry || task.prePlot.plot2Dgeometry) && task.prePlot.abortAfterPlotting
    return
end

%% Build stiffness matrix
t_start = tic;
f = task.f;
if task.useROM && strcmp(task.scatteringCase,'Sweep')
    U_sweep = cell(1,numel(f));
end
if ~(strcmp(task.method,'RT') || strcmp(task.method,'KDT'))
    for i_f = 1:numel(f)
        f_i = f(i_f);
        omega = 2*pi*f_i;
        for m = 1:noDomains
            varCol{m}.omega = omega;
            switch varCol{m}.media
                case 'fluid'
                    varCol{m}.f = f_i;
                    varCol{m}.k = omega/varCol{m}.c_f;
                    varCol{m}.lambda = 2*pi/varCol{m}.k;
                case 'solid'
                    varCol{m}.k = NaN;
            end
        end
        t_freq = tic;
        varCol = getAnalyticSolutions(varCol);
        switch task.method
            case {'IE','ABC'}
                varCol{1}.timeBuildSystem = 0;
                if strcmp(varCol{1}.coreMethod,'SEM')
                    varCol{1} = buildSEMMatrices(varCol{1});
                else
                    for i = 1:noDomains
                        tic          
                        if ~runTasksInParallel
                            fprintf(['\n%-' num2str(stringShift) 's'], ['Building matrices for domain ' num2str(i) ' ... '])
                        end
                        varCol{i} = buildMatrices(varCol{i});
                        varCol{1}.timeBuildSystem = varCol{1}.timeBuildSystem + toc;
                        if ~runTasksInParallel
                            fprintf('using %12f seconds.', toc)
                        end
                    end
                end

                if strcmp(task.method,'IE') && ~strcmp(varCol{1}.coreMethod,'SEM')
                    % Add contribution from infinite elements
                    tic
                    if ~runTasksInParallel
                        fprintf(['\n%-' num2str(stringShift) 's'], 'Building infinite element matrix ... ')
                    end
                    varCol{1} = buildIEmatrix(varCol{1});
                    varCol{1}.timeBuildSystem = varCol{1}.timeBuildSystem + toc;
                    if ~runTasksInParallel
                        fprintf('using %12f seconds.', toc)
                    end
                elseif strcmp(task.method,'ABC')
                    % Add contribution from infinite elements
                    tic
                    if ~runTasksInParallel
                        fprintf(['\n%-' num2str(stringShift) 's'], 'Building ABC matrix ... ')
                    end
                    varCol{1} = addABC(varCol{1}); 

                    varCol{1}.timeBuildSystem = varCol{1}.timeBuildSystem + toc;
                    if ~runTasksInParallel
                        fprintf('using %12f seconds.', toc)
                    end
                end  

                % Apply coupling conditions  
                if noDomains > 1 
                    varCol{2}.p_inc = varCol{1}.p_inc;
                    varCol{2}.dp_inc = varCol{1}.dp_inc;
                    tic
                    if ~runTasksInParallel
                        fprintf(['\n%-' num2str(stringShift) 's'], 'Building coupling matrix ... ')
                    end  
                    for i = 2:noDomains
                        varCol{i} = applyCouplingConditionPatches(varCol{i},varCol{i-1});
                    end
                    if ~runTasksInParallel
                        fprintf('using %12f seconds.', toc)
                    end
                    varCol{1}.timeBuildSystem = varCol{1}.timeBuildSystem + toc;
                end  

                % Apply Neumann conditions
                tic
                if ~runTasksInParallel
                    fprintf(['\n%-' num2str(stringShift) 's'], 'Building RHS vector ... ')
                end
                if task.useROM && noDomains > 1
                    varCol{2}.p_inc_ROM = varCol{1}.p_inc_ROM;
                    varCol{2}.dp_inc_ROM = varCol{1}.dp_inc_ROM;
                end
                for i = 1:min(noDomains,2)
                    varCol{i} = applyNeumannCondition(varCol{i},i==2);
                end
                if ~runTasksInParallel
                    fprintf('using %12f seconds.', toc)
                end
                varCol{1}.timeBuildSystem = varCol{1}.timeBuildSystem + toc;
                
            case 'IENSG'
                tic
                if ~runTasksInParallel
                    fprintf(['\n%-' num2str(stringShift) 's'], 'Building infinite element matrix ... ')
                end
                chimax = varCol{1}.chimax;
                chimin = varCol{1}.chimin;
                if abs(chimax - chimin)/abs(chimax) < 100*eps
                    varCol{1}.r_a = mean([chimax,chimin]);
                    varCol{1} = buildIEmatrix(varCol{1});
                else
                    varCol{1} = infElementsNonSepGeom(varCol{1});  
                end
                varCol{1}.A_K = varCol{1}.Ainf;
                varCol{1} = rmfield(varCol{1},'Ainf');
                varCol{1}.timeBuildSystem = toc;
                if ~runTasksInParallel
                    fprintf('using %12f seconds.', toc)
                end

                % Apply Neumann conditions
                tic
                if ~runTasksInParallel
                    fprintf(['\n%-' num2str(stringShift) 's'], 'Building right hand side vector ... ')
                end

                varCol{1} = applyNeumannCondition(varCol{1});

                varCol{1}.timeBuildSystem = varCol{1}.timeBuildSystem + toc;
                if ~runTasksInParallel
                    fprintf('using %12f seconds.', toc)
                end
            case 'BEM'
                tic
                if ~runTasksInParallel
                    fprintf(['\n%-' num2str(stringShift) 's'], 'Building BEM matrix ... ')
                end
                switch task.formulation(1)
                    case 'G' % Galerkin
                        varCol{1} = buildGBEMmatrix(varCol{1},noDomains > 1);  
                    case 'C' % Collocation
                        varCol{1} = buildCBEMmatrix(varCol{1});  
                end
                if ~runTasksInParallel
                    fprintf('using %12f seconds.', toc)
                end
                if strcmp(task.formulation(end),'C')
                    tic
                    fprintf(['\n%-' num2str(stringShift) 's'], 'Building CHIEF matrix ... ')
                    varCol{1} = buildCHIEFmatrix(varCol{1});
                    if ~runTasksInParallel
                        fprintf('using %12f seconds.', toc)
                    end
                end
                varCol{1}.timeBuildSystem = toc(t_start);
                if noDomains > 1 
                    tic
                    % Solid matrices
                    if ~runTasksInParallel
                        fprintf(['\n%-' num2str(stringShift) 's'], 'Building solid matrix ... ')
                    end
                    varCol{2} = buildMatrices(varCol{2});
                    
                    varCol{1}.timeBuildSystem = varCol{1}.timeBuildSystem + toc;
                    if ~runTasksInParallel
                        fprintf('using %12f seconds.', toc)
                    end
                end
                if noDomains > 2
                    if ~runTasksInParallel
                        fprintf(['\n%-' num2str(stringShift) 's'], 'Building inner fluid matrix ... ')
                    end
                    
                    switch formulation(1)
                        case 'G' % Galerkin
                            varCol{3} = buildGBEMmatrix(varCol{3},noDomains > 1);  
                        case 'C' % Collocation
                            error('Not implemented')
                    end
                    if ~runTasksInParallel
                        fprintf('using %12f seconds.', toc)
                    end
                end

                % Apply coupling conditions  
                if noDomains > 1 
                    tic
                    if ~runTasksInParallel
                        fprintf(['\n%-' num2str(stringShift) 's'], 'Building coupling matrix ... ')
                    end  
                    for i = 2:noDomains
                        if strcmp(varCol{i}.media,'fluid')
                            varCol{i} = applyCouplingCondition_FEBE(varCol{i});
                        end
                    end
                    if ~runTasksInParallel
                        fprintf('using %12f seconds.', toc)
                    end
                    varCol{1}.timeBuildSystem = varCol{1}.timeBuildSystem + toc;
                end  
            case 'BA'
                varCol{1}.timeBuildSystem = 0;
                for i = 1:noDomains
                    tic          
                    if ~runTasksInParallel
                        fprintf(['\n%-' num2str(stringShift) 's'], ['Building matrices for domain ' num2str(i) ' ... '])
                    end
                    varCol{i} = buildMatrices(varCol{i});
                    varCol{i} = buildBAmatrix(varCol{i},i);
                    varCol{1}.timeBuildSystem = varCol{1}.timeBuildSystem + toc;
                    if ~runTasksInParallel
                        fprintf('using %12f seconds.', toc)
                    end
                end
%                 if noDomains > 1 
%                     warning('BA is not implemented in such a way that the displacement and pressure conditions at the interfaces is satisfied')
%                 end
                varCol{1}.timeBuildSystem = toc;
            case 'MFS'
                tic
                if ~runTasksInParallel
                    fprintf(['\n%-' num2str(stringShift) 's'], 'Building MFS matrix ... ')
                end
                varCol{1} = buildMFSmatrix(varCol{1});  
                varCol{1}.timeBuildSystem = toc;
                if ~runTasksInParallel
                    fprintf('using %12f seconds.', varCol{1}.timeBuildSystem)
                end
        end
        switch task.method
            case {'IE','IENSG','BEM','BA','ABC','MFS'}
                [varCol,FF,A0,A1,A2,A4] = collectMatrices(varCol,task);
                if strcmp(task.method,'BA')
                    A = A2;
                else
                    A = A0 + omega*A1 + omega^2*A2 + omega^4*A4;
                end
                
                if ~runTasksInParallel
                    fprintf(['\n%-' num2str(stringShift) 's'], 'Total time building system ... ')
                    fprintf(' was  %12f seconds.', varCol{1}.timeBuildSystem)
                end
                
                %% SOLVE SYSTEM
                tic
                if ~runTasksInParallel
                    fprintf(['\n%-' num2str(stringShift) 's'], ['Creating preconditioner (' task.preconditioner ') ... '])
                end
                switch task.preconditioner
                    case 'ilu'
    %                     [L_A,U_A] = ilu(A,struct('type','ilutp', 'droptol', 1e-4));
            %             [L_A,U_A] = ilu(A,struct('type','ilutp', 'droptol', 1e-3));
                        [L_A,U_A] = ilu(A,struct('type','nofill'));
                    %     [L_A,U_A] = lu(A);
                    case 'SSOR'
                        D_SSOR = spdiags(spdiags(A,0),0,size(A,1),size(A,2));
                        D_SSORinv = spdiags(1./spdiags(A,0),0,size(A,1),size(A,2));
                        F_SSOR = -triu(A);
                        E_SSOR = -tril(A);
                        omega_SSOR = 1.5;
                        L_A = (D_SSOR-omega_SSOR*E_SSOR)*D_SSORinv;
                        U_A = D_SSOR-omega_SSOR*F_SSOR;
%                         P = L_A*U_A;
                    case 'diag'
                        Pinv = spdiags(1./diag(A),0,size(A,1),size(A,2));
                end
                if ~runTasksInParallel
                    fprintf('using %12f seconds.', toc)
                end
                tic
                if ~runTasksInParallel
                    fprintf(['\n%-' num2str(stringShift) 's'], 'Solving system of equations ... ')
                end
                switch task.solver
                    case 'gmres'
                        noRestarts = []; % change this to save memory
                        % noRestarts = 2;
                        [UU,~,~,it1,rv1] = gmres(A,FF,noRestarts,1e-20,1000,L_A,U_A);
                    case 'LU'
                        if ~strcmp(task.preconditioner,'diag')
                            error('not implemented')
                        end
                        if task.useROM
                            dAdomega = A1 + 2*omega*A2 + 4*omega^3*A4;
                            d2Adomega2 = 2*A2 + 12*omega^2*A4;
                            if task.useDGP
                                varCol{1}.A0 = A0;
                                varCol{1}.A1 = A1;
                                varCol{1}.A2 = A2;
                            end
                            UU = zeros(size(FF));
                            dA = decomposition(A*Pinv,'lu');
                            fprintf('using %12f seconds.', toc)
                            fprintf(['\n%-' num2str(stringShift) 's'], 'Computing ROM solution ... ')
                            for i = 1:varCol{1}.noRHSs
                                j = i-1;
                                b = FF(:,i);
                                if j > 0
                                    b = b - j*dAdomega*UU(:,i-1);
                                end
                                if j > 1
                                    b = b - j*(j-1)/2*d2Adomega2*UU(:,i-2);
                                end
                                UU(:,i) = Pinv*(dA\b);
                            end
                        else
                            if i_f == 1 && strcmp(task.formulation,'Sweep')
                                UU = zeros(size(A,1),numel(f));
                            end
                            if strcmp(task.scatteringCase,'Sweep')
                                UU(:,i_f) = Pinv*((A*Pinv)\FF);
                            else
                                UU = Pinv*((A*Pinv)\FF);
                            end
                        end
                    otherwise
                        eval(['[UU,~,~,it1,rv1] = cgs' task.solver '(A,FF,1e-20,1000,L_A,U_A);'])
                end
                dofs = size(A,1);
                varCol{1}.dofs = dofs;
                if task.clearGlobalMatrices
                    varCol = rmfields(varCol,{'A_K','A_M','A_C','FF','Ainf','Ainf1','Ainf2'});
                    clear A FF A0 A1 A2 A4 Pinv
                end
                if ~runTasksInParallel
                    fprintf('using %12f seconds.', toc)
                end
                varCol{1}.timeSolveSystem = toc;
                if task.computeCondNumber && (size(A,1) == size(A,2))
                    fprintf(['\n%-' num2str(stringShift) 's'], 'Calculating condition number ... ')
                    rng('default') % for reproducibility in condest
                    condNumber = condest(A);
        %             condNumber = cond(full(A))
                    fprintf('using %12f seconds.', toc)
                    fprintf('\nCondition number = %d', condNumber)
                else
                    condNumber = NaN;
                end
                task.results.cond_number = condNumber;

                if ~runTasksInParallel
                    fprintf('\nNumber of degrees of freedom = %d', dofs)
                end
                % fprintf('\nMemory ratio = %f', ((fluid.degree(1)+1)^6*varCol{1}.noElems)/nnz(A_fluid_o))
        end
        if task.useROM && strcmp(task.scatteringCase,'Sweep')
            [varCol,U_sweep{i_f}] = postProcessSolution(varCol,UU);
        end
        if strcmp(task.scatteringCase,'Sweep') && numel(f) > 1
            fprintf('\nTotal time spent on frequency %d of %d: %12f\n', i_f, numel(f), toc(t_freq))  
        end
    end
    if numel(f) > 1 && ~task.useROM && (task.calculateSurfaceError || task.calculateVolumeError)
        tic
        fprintf(['\n%-' num2str(stringShift) 's'], 'Calculating surface error ... ')
        progressBars = numel(f) > 1 && task.progressBars;
        nProgressStepSize = ceil(numel(f)/10);
        if progressBars
            ppm = ParforProgMon('Calculating surface error: ', numel(f));
        else
            ppm = NaN;
        end
    else
        progressBars = false;
        ppm = NaN;
    end

    if ~task.useROM && (task.calculateSurfaceError || task.calculateVolumeError)
        for i_f = 1:numel(f)
            if progressBars
                ppm.increment();
            end
            f_i = f(i_f);
            omega = 2*pi*f_i;
            for m = 1:noDomains
                varCol{m}.omega = omega;
                switch varCol{m}.media
                    case 'fluid'
                        varCol{m}.f = f_i;
                        varCol{m}.k = omega/varCol{m}.c_f;
                        varCol{m}.lambda = 2*pi/varCol{m}.k;
                    case 'solid'
                        varCol{m}.k = NaN;
                end
            end
            varCol = getAnalyticSolutions(varCol);
            if i_f == 1
                task.results.energyError = zeros(1,size(f,2));
                task.results.L2Error = zeros(1,size(f,2));
                task.results.H1Error = zeros(1,size(f,2));
                task.results.H1sError = zeros(1,size(f,2));
                task.results.surfaceError = zeros(1,size(f,2));
            end
            varCol = postProcessSolution(varCol,UU(:,i_f));
            printLog = numel(f) == 1 && ~runTasksInParallel;
            printLog = true && ~runTasksInParallel;
            [L2Error, H1Error, H1sError, energyError, surfaceError] ...
                = calculateErrors(task, varCol, printLog, stringShift);
            task.results.surfaceError(i_f) = surfaceError;
            task.results.energyError(i_f) = energyError;
            task.results.L2Error(i_f) = L2Error;
            task.results.H1Error(i_f) = H1Error;
            task.results.H1sError(i_f) = H1sError;
        end
    end
    if numel(f) > 1 && ~task.useROM && (task.calculateSurfaceError || task.calculateVolumeError)
        fprintf('using %12f seconds.', toc)   
    end
else
    varCol{1}.U = [];
end  

%% Compute scattered pressure   
if task.calculateFarFieldPattern && ~task.useROM
    omega = 2*pi*f;
    varCol{1}.omega = omega;
    varCol{1}.k = omega/varCol{1}.c_f;
    varCol{1}.lambda = 2*pi./varCol{1}.k;
    varCol{1}.f = f;
    varCol = getAnalyticSolutions(varCol);
    varCol = postProcessSolution(varCol,UU);
    task = calculateTS(varCol,task,runTasksInParallel,stringShift);
end
if ~task.useROM
    varCol{1}.k = (2*pi*f)/varCol{1}.c_f;
    varCol{1}.f = f;
    varCol{1}.omega = f/(2*pi);
end
if task.clearGlobalMatrices
    clear UU U
end

%% Calculate errors (if analyticSolutionExist) and plot result in Paraview
if ~task.useROM && ~strcmp(task.method,'RT')
    switch task.scatteringCase
        case {'BI', 'Sweep','Ray'}    
            tic
            if task.plotResidualError && varCol{1}.analyticSolutionExist && strcmp(task.formulation,'GCBIE')
%                 close all
                fprintf(['\n%-' num2str(stringShift) 's'], 'Plotting zeros of residual ... ')
                plotGalerkinResidual(varCol{1});
                fprintf('using %12f seconds.', toc)
                axis equal
                axis off
    %             set(gca, 'Color', 'none');
%                 title('Error in Galerkin solution')
%                 colorbar
%                 caxis([-7,0])
%                 camlight
                colorbar off
                savefig([resultsFolder '/' task.saveName])
    %             n_xi = fluid.number(1);
    %             n_eta = fluid.number(2);
    %             [cg_xi, grev_xi] = CauchyGalerkin(fluid.degree(1), n_xi, fluid.knots{1});
    %             [cg_eta, grev_eta] = CauchyGalerkin(fluid.degree(2), n_eta, fluid.knots{2});
    %             % keyboard
    %             n_cp = varCol{1}.noDofs - length(dofsToRemove);
    %             
    %             cp_cg = zeros(n_cp,3);
    %             cp_grev = zeros(n_cp,3);
    %             counter = 1;
    %             counter2 = 1;
    %             for j = 1:n_eta
    %                 eta_cg = cg_eta(j);
    %                 eta_grev = grev_eta(j);
    %                 for i = 1:n_xi
    %                     if ~any(dofsToRemove == counter)
    %                         xi_cg = cg_xi(i);        
    %                         xi_grev = grev_xi(i);            
    %                         cp_cg(counter2,:) = evaluateNURBS(fluid, [xi_cg, eta_cg]); 
    %                         cp_grev(counter2,:) = evaluateNURBS(fluid, [xi_grev, eta_grev]); 
    %                         counter2 = counter2 + 1;
    %                     end
    %                     counter = counter + 1;
    %                 end
    %             end
    %             hold on
    %             h1 = plot3(cp_cg(:,1),cp_cg(:,2),cp_cg(:,3), '*','color','red');
    %             h2 = plot3(cp_grev(:,1),cp_grev(:,2),cp_grev(:,3), '*','color','blue');
    %             legend([h1, h2], {'CG','Grev'})
    %             hold off
    %             savefig([resultsFolderName '/' task.saveName '_surfPlot_mesh' num2str(M) '_formulation_' formulation '_degree' num2str(max(fluid.degree)) '.fig'])
            end
            para = task.para;
            if para.plotResultsInParaview
                tic
                if ~runTasksInParallel
                    fprintf(['\n%-' num2str(stringShift) 's'], 'Post-processing ... ')
                end
                resultsFolderNameParaview = [resultsFolder '/paraviewResults'];
                if ~exist(resultsFolderNameParaview, 'dir')
                    mkdir(resultsFolderNameParaview);
                end          
                vtfFileName = [resultsFolderNameParaview '/' task.saveName];
                if isempty(para.name)
                    para.name = vtfFileName;
                end


                para.plotArtificialBoundary = para.plotArtificialBoundary && (strcmp(task.method,'IE') || strcmp(task.method,'ABC') || strcmp(task.method,'PML'));
                M = task.M;
                para.extraXiPts = eval(para.extraXiPts);
                para.extraEtaPts = eval(para.extraEtaPts);
                para.extraZetaPts = eval(para.extraZetaPts);
                testFun = @(v) -analytic(v) - P_inc(v);
                varCol{1}.testFun = testFun;
                if strcmp(task.method, 'KDT')
                    varCol{1}.U = zeros(varCol{1}.noDofs,1);
                end
        
                createParaviewFiles(varCol, 'para_options', para);

                % plot artificial boundary
%                 if para.plotArtificialBoundary
%                     plotModelInParaview(subNURBS(varCol{1}.nurbs,'at',[0,0;0,0;0,1]), para, 0, 'artificialBoundary')
%                     if strcmp(BC,'SHBC')
%                         createParaviewFiles(subNURBS(varCol{1}.nurbs,'at',[0,0;0,0;1,0]), para, 1, 'scatterer')
%                     end
%                 end
            end
    end
end
varCol{1}.tot_time = toc(t_start);

fprintf('\nTotal time spent on task: %12f', varCol{1}.tot_time)  

if task.storeFullVarCol
    varColTemp = varCol;
    if task.useROM
        varColTemp{1}.U_sweep = U_sweep;
        if task.useDGP
            varColTemp{1}.A0 = A0;
            varColTemp{1}.A1 = A1;
            varColTemp{1}.A2 = A2;
        end
    end
    task.varCol = varColTemp;
else
    task.varCol = extractVarColFields(task,varCol);
end






