function task = main_sub(task,loopParameters,printLog,resultsFolder)

stringShift = 40;
task = controlTask(task);
task.saveName = defineFolderAndFilenames(task,loopParameters);
task.resultsFolder = resultsFolder;
if printLog
    fprintf('\n\nRunning the case: %s', task.saveName)
end

task = setAndDefineParameters(task);

if printLog
    fprintf(['\n%-' num2str(stringShift) 's'], 'Extracting CAD data ... ')
    tic
end
task = createNURBSmesh(task);
task = collectVariables(task);
if printLog
    fprintf('using %12f seconds.', toc)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot geometry
% keyboard
if task.prePlot.plot3Dgeometry || task.prePlot.plot2Dgeometry
    tic
    fprintf(['\n%-' num2str(stringShift) 's'], 'Plotting geometry ... ')
    task = plotMeshAndGeometry(task);
    fprintf('using %12f seconds.', toc)
end

%% Build connectivity
tic
if printLog
    fprintf(['\n%-' num2str(stringShift) 's'], 'Generating IGA mesh ... ')
end
for i = 1:task.noDomains
    task.varCol{i} = generateIGAmesh(task.varCol{i});
end

%% Find nodes to remove due to gluing and homogeneous dirichlet conditions
for i = 1:task.noDomains
    task.varCol{i} = findDofsToRemove(task.varCol{i});
end
if printLog
    fprintf('using %12f seconds.', toc)
end

%% Compute derived quantities
tic
if printLog
    fprintf(['\n%-' num2str(stringShift) 's'], 'Computing derived quantities ... ')
end
task = computeDerivedQuantities(task);

%% Check NURBS compatibility
if task.misc.checkNURBSweightsCompatibility
    checkNURBSweightsCompatibility(task);
end
if ~isempty(task.prePlot.QoI)
    task.results.QoIError = integrateFunc(task, task.prePlot.QoI, task.prePlot.QoI_ref);
    fprintf('\nQoI error is %.15g\n',task.results.QoIError)
end

if printLog
    fprintf('using %12f seconds.', toc)
    fprintf('\nTotal number of elements = %d', task.totNoElems)
    fprintf('\nFinite element dofs = %d', task.FEdofs)
    fprintf('\nNumber of elements per wavelength = %.2g', min(task.varCol{1}.nepw(:)))
end

if (task.prePlot.plot3Dgeometry || task.prePlot.plot2Dgeometry || task.misc.preProcessOnly) && task.prePlot.abortAfterPlotting
    return
end

%% Build stiffness matrix
t_start = tic;
omega = task.misc.omega;
if task.rom.useROM && strcmp(task.misc.scatteringCase,'Sweep')
    task.U_sweep = cell(1,numel(omega));
end
if ~(strcmp(task.misc.method,'RT') || strcmp(task.misc.method,'KDT'))
    for i_o = 1:numel(omega)
        omega_i = omega(i_o);
        t_freq = tic;
        task.misc.omega = omega_i;
        task = getAnalyticSolutions(task);
        switch task.misc.method
            case {'IE','ABC','IENSG','PML'}
                task.timeBuildSystem = 0;
                for i_domain = 1:task.noDomains
                    if ~(strcmp(task.misc.method,'IENSG') && task.iem.boundaryMethod && i_domain == 1)
                        tic          
                        if printLog
                            fprintf(['\n%-' num2str(stringShift) 's'], ['Building matrices for domain ' num2str(i_domain) ' ... '])
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
                        if printLog
                            fprintf('using %12f seconds.', toc)
                        end
                    end
                end
                if ~strcmp(task.misc.coreMethod,'SEM')
                    if strcmp(task.misc.method,'IE')
                        % Add contribution from infinite elements
                        tic
                        if printLog
                            fprintf(['\n%-' num2str(stringShift) 's'], 'Building infinite element matrix ... ')
                        end
                        task = buildIEmatrix(task);
                        
                        task.timeBuildSystem = task.timeBuildSystem + toc;
                        if printLog
                            fprintf('using %12f seconds.', toc)
                        end
                    elseif strcmp(task.misc.method,'IENSG')
                        tic
                        if printLog
                            fprintf(['\n%-' num2str(stringShift) 's'], 'Building infinite element matrix ... ')
                        end
                        chimax = task.varCol{1}.chimax;
                        chimin = task.varCol{1}.chimin;
                        if abs(chimax - chimin)/abs(chimax) < 100*eps % Exploit tensor product structure of IEM
                            task.misc.r_a = mean([chimax,chimin]);
                            task = buildIEmatrix(task);
                        else
                            task = infElementsNonSepGeom(task);  
                        end
                        if printLog
                            fprintf('using %12f seconds.', toc)
                        end
                    elseif strcmp(task.misc.method,'ABC')
                        % Add contribution from infinite elements
                        tic
                        if printLog
                            fprintf(['\n%-' num2str(stringShift) 's'], 'Building ABC matrix ... ')
                        end
                        task.varCol{1} = addABC(task.varCol{1}); 

                        task.timeBuildSystem = task.timeBuildSystem + toc;
                        if printLog
                            fprintf('using %12f seconds.', toc)
                        end
                    end  
                    % Apply coupling conditions  
                    if task.noDomains > 1 
                        tic
                        if printLog
                            fprintf(['\n%-' num2str(stringShift) 's'], 'Building coupling matrix ... ')
                        end
                        for i_domain = 2:task.noDomains
                            task = applyCouplingConditionPatches(task,i_domain);
                        end
                        if printLog
                            fprintf('using %12f seconds.', toc)
                        end
                        task.timeBuildSystem = task.timeBuildSystem + toc;
                    end

                    % Apply Neumann conditions
                    tic
                    if printLog
                        fprintf(['\n%-' num2str(stringShift) 's'], 'Building RHS vector ... ')
                    end
                    for i_domain = 1:min(task.noDomains,2)
                        task = applyNeumannCondition(task,i_domain);
                        if task.misc.symmetric && ~task.rom.useROM && strcmp(task.varCol{i_domain}.media,'solid')
                            task.varCol{i_domain}.FF = omega_i^2*task.varCol{i_domain}.FF;
                        end
                    end
                    if printLog
                        fprintf('using %12f seconds.', toc)
                    end
                    task.timeBuildSystem = task.timeBuildSystem + toc;
                end
            case 'BEM'
                tic
                if printLog
                    fprintf(['\n%-' num2str(stringShift) 's'], 'Building BEM matrix ... ')
                end
                switch task.misc.formulation(1)
                    case 'G' % Galerkin
                        task = buildGBEMmatrix(task,1);  
                    case 'C' % Collocation
                        task = buildCBEMmatrix(task);  
                end
                if printLog
                    fprintf('using %12f seconds.', toc)
                end
                if strcmp(task.misc.formulation(end),'C')
                    tic
                    if printLog
                        fprintf(['\n%-' num2str(stringShift) 's'], 'Building CHIEF matrix ... ')
                    end
                    task = buildCHIEFmatrix(task);
                    if printLog
                        fprintf('using %12f seconds.', toc)
                    end
                end
                task.timeBuildSystem = toc(t_start);
                if task.noDomains > 1 
                    tic
                    % Solid matrices
                    if printLog
                        fprintf(['\n%-' num2str(stringShift) 's'], 'Building solid matrix ... ')
                    end
                    task = buildMatrices(task,2);
                    
                    task.timeBuildSystem = task.timeBuildSystem + toc;
                    if printLog
                        fprintf('using %12f seconds.', toc)
                    end
                end
                if task.noDomains > 2
                    if printLog
                        fprintf(['\n%-' num2str(stringShift) 's'], 'Building inner fluid matrix ... ')
                    end
                    
                    switch formulation(1)
                        case 'G' % Galerkin
                            task = buildGBEMmatrix(task,3);  
                        case 'C' % Collocation
                            error('Not implemented')
                    end
                    if printLog
                        fprintf('using %12f seconds.', toc)
                    end
                end

                % Apply coupling conditions  
                if task.noDomains > 1 
                    tic
                    if printLog
                        fprintf(['\n%-' num2str(stringShift) 's'], 'Building coupling matrix ... ')
                    end  
                    for i = 2:task.noDomains
                        if strcmp(task.varCol{i}.media,'fluid')
                            task.varCol{i} = applyCouplingCondition_FEBE(task.varCol{i});
                        end
                    end
                    if printLog
                        fprintf('using %12f seconds.', toc)
                    end
                    task.timeBuildSystem = task.timeBuildSystem + toc;
                end  
            case 'BA'
                task.timeBuildSystem = 0;
                for i = 1:task.noDomains
                    tic          
                    if printLog
                        fprintf(['\n%-' num2str(stringShift) 's'], ['Building matrices for domain ' num2str(i) ' ... '])
                    end
                    task = buildMatrices(task,i);
                    task = buildBAmatrix(task,i);
                    task.timeBuildSystem = task.timeBuildSystem + toc;
                    if printLog
                        fprintf('using %12f seconds.', toc)
                    end
                end
%                 if task.noDomains > 1 
%                     warning('BA is not implemented in such a way that the displacement and pressure conditions at the interfaces is satisfied')
%                 end
                task.timeBuildSystem = toc;
            case 'MFS'
                tic
                if printLog
                    fprintf(['\n%-' num2str(stringShift) 's'], 'Building MFS matrix ... ')
                end
                task = buildMFSmatrix(task);  
                task.timeBuildSystem = toc;
                if printLog
                    fprintf('using %12f seconds.', task.timeBuildSystem)
                end
        end
        switch task.misc.method
            case {'IE','IENSG','BEM','BA','ABC','MFS','PML'}
                tic
                if printLog
                    fprintf(['\n%-' num2str(stringShift) 's'], 'Collecting matrices ... ')
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
                if printLog
                    fprintf('using %12f seconds.', task.timeCollectingMatrices)
                end
                
                if printLog
                    fprintf(['\n%-' num2str(stringShift) 's'], 'Total time building system ... ')
                    fprintf(' was  %12f seconds.', task.timeBuildSystem+task.timeCollectingMatrices)
                end
                
                %% SOLVE SYSTEM
                tic
                if printLog
                    fprintf(['\n%-' num2str(stringShift) 's'], ['Creating preconditioner (' task.sol.preconditioner ') ... '])
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
                if printLog
                    fprintf('using %12f seconds.', toc)
                end
                tic
                if printLog
                    fprintf(['\n%-' num2str(stringShift) 's'], 'Solving system of equations ... ')
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
                            fprintf(['\n%-' num2str(stringShift) 's'], 'Computing ROM solution ... ')
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
                if printLog
                    fprintf('using %12f seconds.', toc)
                end
                task.timeSolveSystem = toc;
                if task.misc.computeCondNumber && (size(A,1) == size(A,2))
                    if printLog
                        fprintf(['\n%-' num2str(stringShift) 's'], 'Calculating condition number ... ')
                    end
                    rng('default') % for reproducibility in condest
                    condNumber = condest(A);
        %             condNumber = cond(full(A))
                    if printLog
                        fprintf('using %12f seconds.', toc)
                        fprintf('\nCondition number = %d', condNumber)
                    end
                else
                    condNumber = NaN;
                end
                task.results.cond_number = condNumber;

                if printLog
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
        if printLog && strcmp(task.misc.scatteringCase,'Sweep') && numel(omega) > 1
            fprintf('\nTotal time spent on frequency %d of %d: %12f\n', i_o, numel(omega), toc(t_freq))  
        end
    end
    if numel(omega) > 1 && ~task.rom.useROM && (task.err.calculateSurfaceError || task.err.calculateVolumeError)
        tic
        if printLog
            fprintf(['\n%-' num2str(stringShift) 's'], 'Calculating surface error ... ')
        end
        progressBars = numel(omega) > 1 && task.misc.progressBars;
        nProgressStepSize = ceil(numel(omega)/10);
        if progressBars
            try
                ppm = ParforProgMon('Calculating surface error: ', numel(omega));
            catch
                progressBars = false;
                ppm = NaN;
            end
        else
            ppm = NaN;
        end
    else
        progressBars = false;
        ppm = NaN;
    end

    if ~task.rom.useROM
        task = postProcessSolution(task);
    end
    if ~task.rom.useROM && (task.err.calculateSurfaceError || task.err.calculateVolumeError)
        task.results.energyError = zeros(1,numel(omega));
        task.results.L2Error = zeros(1,numel(omega));
        task.results.H1Error = zeros(1,numel(omega));
        task.results.H1sError = zeros(1,numel(omega));
        task.results.surfaceError = zeros(1,numel(omega));
        for i_o = 1:numel(omega)
            if progressBars
                ppm.increment();
            end
            task.misc.omega = omega(i_o);
            task = getAnalyticSolutions(task);

%             printLog = numel(omega) == 1 && printLog;
            [L2Error, H1Error, H1sError, energyError, surfaceError] = calculateErrors(task, printLog, stringShift);
            task.results.surfaceError(i_o) = surfaceError;
            task.results.energyError(i_o) = energyError;
            task.results.L2Error(i_o) = L2Error;
            task.results.H1Error(i_o) = H1Error;
            task.results.H1sError(i_o) = H1sError;
        end
    end
    if numel(omega) > 1 && ~task.rom.useROM && (task.err.calculateSurfaceError || task.err.calculateVolumeError)
        fprintf('using %12f seconds.', toc)   
    end
else
    task.varCol{1}.U = [];
end  
task.misc.omega = omega;

%% Compute scattered pressure   
if task.ffp.calculateFarFieldPattern && ~task.rom.useROM
    task = getAnalyticSolutions(task);
    task = calculateTS(task,printLog,stringShift);
    task = computeDerivedFFPquantities(task,task.ffp.p_h);
end
%% Calculate errors (if analyticSolutionExist) and plot result in Paraview
if ~task.rom.useROM && ~strcmp(task.misc.method,'RT')
    tic
    if task.misc.plotResidualError && task.analyticSolutionExist && strcmp(task.misc.formulation,'GCBIE')
%                 close all
        fprintf(['\n%-' num2str(stringShift) 's'], 'Plotting zeros of residual ... ')
        plotGalerkinResidual(task.varCol{1});
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
%             n_cp = task.varCol{1}.noDofs - length(dofsToRemove);
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
        task.ffp.alpha_s = task.ffp.alpha_s(task.para.i_MS);
        task = getAnalyticSolutions(task);
        if isempty(para.name)
            resultsFolderNameParaview = [resultsFolder '/paraviewResults'];
            if ~exist(resultsFolderNameParaview, 'dir')
                mkdir(resultsFolderNameParaview);
            end          
            vtfFileName = [resultsFolderNameParaview '/' task.saveName];
            para.name = vtfFileName;
        end

        para.extraXiPts = eval(para.extraXiPts);
        para.extraEtaPts = eval(para.extraEtaPts);
        para.extraZetaPts = eval(para.extraZetaPts);
        if strcmp(task.misc.method, 'KDT')
            task.varCol{1}.U = zeros(task.varCol{1}.noDofs,1);
        end
        if para.plotFullDomain
            tic
            if printLog
                fprintf(['\n%-' num2str(stringShift) 's'], 'Creating volumetric paraview files for all domains ... ')
            end
            createParaviewFiles(task, 'para_options', para);
        end

        for i_v = 1:numel(task.varCol)
            if ~isempty(para.plotSubsets) && ~task.varCol{i_v}.boundaryMethod
                name = para.name;
                for i = 1:numel(para.plotSubsets)
                    bdryName = para.plotSubsets{i};
                    if printLog
                        fprintf(['\n%-' num2str(stringShift) 's'], ['Creating paraview files for ' bdryName ' in domain ' num2str(i_v) ' ... '])
                    end
                    [~,setFound] = findSet(task.varCol{i_v}.geometry.topologysets.set,bdryName);
                    if setFound
                        taskBdry = task;
                        taskBdry.varCol{i_v} = rmfield(taskBdry.varCol{i_v},'geometry');
                        varColBdry = meshBoundary(task.varCol{i_v},bdryName);
                        taskBdry.varCol{i_v}.nurbs = varColBdry.nurbs;
                        taskBdry = collectVariables(taskBdry);
                        taskBdry.varCol{i_v} = findDofsToRemove(generateIGAmesh(taskBdry.varCol{i_v}));
    %                     fieldNames = fields(varColBdry);
    %                     for j = 1:numel(fieldNames)
    %                         taskBdry.varCol{i_v}.(fieldNames{j}) = varColBdry.(fieldNames{j});
    %                     end
                        nodes = varColBdry.nodes;

                        d_f = task.varCol{i_v}.dimension;
                        nodesField = zeros(1,d_f*numel(nodes));
                        for j = 1:d_f
                            nodesField(j:d_f:end) = d_f*(nodes-1)+j;
                        end
                        taskBdry.varCol{i_v}.U = task.varCol{i_v}.U(nodesField,:);

                        para.name = [name, '_', bdryName];

                        createParaviewFiles(taskBdry, 'para_options', para);
                    else
                        switch bdryName
                            case {'xy','yz','xz'}
                                warning(['The set ' bdryName ' was not found for domain ' num2str(i_v) ', and is thus not plotted in paraview. For xy, yz and xz; setting msh.explodeNURBS = true might help.'])
                            otherwise
                                if ~(strcmp(bdryName,'Gamma_a') && i_v > 1)
                                    warning(['The set ' bdryName ' was not found for domain ' num2str(i_v) ', and is thus not plotted in paraview.'])
                                end
                        end
                    end
                end
            end
        end
    end
end
task.varCol{1}.tot_time = toc(t_start);

if printLog
    fprintf('\nTotal time spent on task: %12f', task.varCol{1}.tot_time)  
end

if ~task.misc.storeFullVarCol
    task.varCol = rmfields(task.varCol,getAddedFields());
end






