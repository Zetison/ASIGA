function task = main_sub(task,loopParameters,runTasksInParallel,resultsFolder)

stringShift = 40;
varCol = extractTaskData(task);
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
if task.useROM && strcmp(scatteringCase,'Sweep')
    U_sweep = cell(1,numel(f));
end
f = task.f;
if strcmp(task.method,'RT') || strcmp(task.method,'KDT')
    omega = 2*pi*f;
    varCol{1}.omega = omega;
    varCol{1}.k = omega/varCol{1}.c_f;
    varCol{1}.lambda = 2*pi./varCol{1}.k;
    varCol{1}.f = f;
    varCol = getAnalyticSolutions(varCol);
else
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
        k_i = varCol{1}.k;
        varCol = getAnalyticSolutions(varCol);
        rho = varCol{1}.rho;
        switch task.method
            case {'IE','ABC'}
                tic  
                if ~runTasksInParallel
                    fprintf(['\n%-' num2str(stringShift) 's'], 'Building outer fluid matrix ... ')
                end
                if strcmp(varCol{1}.coreMethod,'SEM')
                    varCol{1} = buildSEMMatrices(varCol{1});
                else
                    varCol{1} = buildMatrices(varCol{1});
                end
                varCol{1}.timeBuildSystem = toc;
                if ~runTasksInParallel
                    fprintf('using %12f seconds.', toc)
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
                    tic        
                    % Inner fluid   
                    if ~runTasksInParallel
                        fprintf(['\n%-' num2str(stringShift) 's'], 'Building inner fluid matrix ... ')
                    end
                    varCol{3} = buildMatrices(varCol{3});

                    shift = varCol{3}.noDofs;
                    if ~runTasksInParallel
                        fprintf('using %12f seconds.', toc)
                    end
                    varCol{1}.timeBuildSystem = varCol{1}.timeBuildSystem + toc;
                end

                % Apply coupling conditions  
                if noDomains > 1 
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
                switch formulation(1)
                    case 'G' % Galerkin
                        varCol{1} = buildGBEMmatrix(varCol{1},noDomains > 1);  
                    case 'C' % Collocation
                        varCol{1} = buildCBEMmatrix(varCol{1});  
                end
                noDofs_tot = varCol{1}.noDofs;
                if ~runTasksInParallel
                    fprintf('using %12f seconds.', toc)
                end
                if strcmp(formulation(end),'C')
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
                    A_coupling = applyCouplingCondition_FEBE(varCol{1});
                    d = 3;
                    solidSurfaceDofs = d*varCol{1}.noDofs;
                    shift2 = shift+varCol{2}.noDofs;
                    A(shift2-solidSurfaceDofs+1:shift2,shift2+1:shift2+varCol{1}.noDofs) = A_coupling.';
                    if ~runTasksInParallel
                        fprintf('using %12f seconds.', toc)
                    end
                    varCol{1}.timeBuildSystem = varCol{1}.timeBuildSystem + toc;
                end

                if noDomains > 2 
                    tic
                    if ~runTasksInParallel
                        fprintf(['\n%-' num2str(stringShift) 's'], 'Building second coupling matrix ... ')
                    end
                    A_coupling = applyCouplingCondition_FEBE(varCol{3});
                    noDofsInner = varCol{3}.noDofsInner;
                    solidSurfaceDofs = d*(varCol{3}.noDofs-noDofsInner);
                    shift2 = varCol{3}.noDofs;
                    A(shift2+1:shift2+solidSurfaceDofs,noDofsInner+1:shift2) = ...
                        -A_coupling(noDofsInner+1:end,d*noDofsInner+1:end).';
                    if ~runTasksInParallel
                        fprintf('using %12f seconds.', toc)
                    end
                    varCol{1}.timeBuildSystem = varCol{1}.timeBuildSystem + toc;
                end  
            case 'BA'
                tic
                if ~runTasksInParallel
                    fprintf(['\n%-' num2str(stringShift) 's'], 'Building outer fluid matrix ... ')
                end
                varCol{1} = buildBAmatrix(varCol{1});
                varCol{1}.timeBuildSystem = toc;
                if ~runTasksInParallel
                    fprintf('using %12f seconds.', varCol{1}.timeBuildSystem)
                end
                if noDomains > 1 
                    error('BA is not implemented in such a way that the displacement and pressure conditions at the interfaces is satisfied')
                end
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
                Aindices = cell(1,noDomains);
                noDofs_tot = 0;
                dofsToRemove = [];
                for i = noDomains:-1:1
                    Aindices{i} = noDofs_tot+(1:varCol{i}.noDofs);
                    dofsToRemove = [dofsToRemove (varCol{i}.dofsToRemove+noDofs_tot)];
                    noDofs_tot = noDofs_tot + varCol{i}.noDofs;
                end
                if strcmp(task.method,'IE') || strcmp(task.method,'IENSG')
                    AindicesInf = noDofs_tot+(1:varCol{1}.noDofs_new);
                    noDofs_tot = noDofs_tot - varCol{1}.noDofs + varCol{1}.noDofs_new;
                end
                % Collect all matrices
                A = sparse(noDofs_tot,noDofs_tot);
                FF = zeros(noDofs_tot,varCol{1}.noRHSs);
                for i = 1:noDomains
%                     j = noDomains-i+1;
                    A(Aindices{i},Aindices{i}) = varCol{i}.A_K - varCol{i}.k*varCol{i}.A_M; 
                    FF(Aindices{i},:) = varCol{i}.FF;
                end
                A(AindicesInf,AindicesInf) = A(AindicesInf,AindicesInf) + varCol{1}.Ainf; 
                A(dofsToRemove,:) = [];
                A(:,dofsToRemove) = [];
                FF(dofsToRemove,:) = [];
                
                
                if ~runTasksInParallel
                    fprintf(['\n%-' num2str(stringShift) 's'], 'Total time building system ... ')
                    fprintf(' was  %12f seconds.', varCol{1}.timeBuildSystem)
                end
                
                %% SOLVE SYSTEM
                if ~strcmp(task.solver,'LU')
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
                    end
                    if ~runTasksInParallel
                        fprintf('using %12f seconds.', toc)
                    end
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
                        
                                
                        if task.useROM
                            A = A_K - k_i^2*A_M + k_i^2*A_2 + k_i*A_1 + A_0;
                            dAdk = -2*k_i*A_M + 2*k_i*A_2 + A_1;
                            d2Adk2 = -2*A_M + 2*A_2;
                            if task.useDGP
                                varCol{1}.A_gamma_a = A_2;
                                varCol{1}.A2_gamma_a = A_1;
                                varCol{1}.A3_gamma_a = A_0;
                            end
                            UU = zeros(size(FF));
                            dA = decomposition(A,'lu');
                            fprintf('using %12f seconds.', toc)
                            fprintf(['\n%-' num2str(stringShift) 's'], 'Computing ROM solution ... ')
                            for i = 1:noVecs
                                j = i-1;
                                b = FF(:,i);
                                if j > 0
                                    b = b - j*dAdk*UU(:,i-1);
                                end
                                if j > 1
                                    b = b - j*(j-1)/2*d2Adk2*UU(:,i-2);
                                end
                                UU(:,i) = dA\b;
                            end
                            fprintf('using %12f seconds.', toc)
                        else
                            A = varCol{1}.A_K - varCol{1}.k^2*varCol{1}.A_M + varCol{1}.Ainf;
                            if strcmp(task.method,'MFS') && strcmp(formulation,'SS')
                                Pinv = diag(1./max(abs(A)));
                                UU = diag(Pinv).*((A*Pinv)\FF);
                            else
                                UU = A\FF;
                            end
                        end
                    otherwise
                        eval(['[UU,~,~,it1,rv1] = cgs' task.solver '(A,FF,1e-20,1000,L_A,U_A);'])
                end
                if ~runTasksInParallel
                    fprintf('using %12f seconds.', toc)
                end
                varCol{1}.timeSolveSystem = toc;
                if computeCondNumber && (size(A,1) == size(A,2))
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

                dofs = size(A,1);
                if ~runTasksInParallel
                    fprintf('\nNumber of degrees of freedom = %d', dofs)
                end
                % fprintf('\nMemory ratio = %f', ((fluid.degree(1)+1)^6*varCol{1}.noElems)/nnz(A_fluid_o))
                varCol{1}.dofs = dofs;
                [varCol,U,Uc] = postProcessSolution(varCol,task,UU);
                if varCol{1}.boundaryMethod
                    varCol{1}.tau = computeTau(varCol{1});
                end
        end
        if task.useROM && strcmp(scatteringCase,'Sweep')
            U_sweep{i_f} = Uc{1}(1:varCol{1}.noDofs,:);
            U_sweep2{i_f} = Uc{1};
            if strcmp(BC,'SSBC') || strcmp(BC,'NNBC')
                error('not implemented due to noDofs')
            end
        end
        if ~task.useROM
            if i_f == 1
                task.results.energyError = zeros(1,size(f,2));
                task.results.L2Error = zeros(1,size(f,2));
                task.results.H1Error = zeros(1,size(f,2));
                task.results.H1sError = zeros(1,size(f,2));
                task.results.surfaceError = zeros(1,size(f,2));
            end
            [L2Error, H1Error, H1sError, energyError, surfaceError] ...
                = calculateErrors(task, varCol, Uc, runTasksInParallel, stringShift, i_f);
            task.results.surfaceError(i_f) = surfaceError;
            task.results.energyError(i_f) = energyError;
            task.results.L2Error(i_f) = L2Error;
            task.results.H1Error(i_f) = H1Error;
            task.results.H1sError(i_f) = H1sError;
        end
        if strcmp(scatteringCase,'Sweep') && size(k_i,2) > 1
            fprintf('\nTotal time spent on frequency: %12f\n', toc(t_freq))  
        end
    end
end
if ~task.useROM
    varCol{1}.k = (2*pi*f)/varCol{1}.c_f;
end

%% Compute scattered pressure    
tic
if calculateFarFieldPattern && ~task.useROM
    if ~runTasksInParallel
        fprintf(['\n%-' num2str(stringShift) 's'], 'Computing far-field pattern ... ')
    end
    v = getFarFieldPoints(task.alpha,task.beta,task.r);

    switch task.method
        case {'IE','ABC','IENSG','BA','BEM'}
            p_h = calculateScatteredPressure(varCol, Uc, v, 0, plotFarField);
        case 'MFS'
            p_h = calculateScatteredPressureMFS(varCol{1}, Uc{1}, v, plotFarField);
        case 'KDT'
            switch varCol{1}.coreMethod
                case 'linear_FEM'
                    noElems = varCol{1}.noElems;
                    element = varCol{1}.element;
                    tri = NaN(size(element,1),2,3);
                    P = varCol{1}.controlPts;
                    Eps = 1e2*eps;
                    for e = 1:noElems
                        sctr = element(e,:);
                        P1 = P(sctr(1),:);
                        P2 = P(sctr(2),:);
                        P3 = P(sctr(3),:);
                        P4 = P(sctr(4),:);
                        tri_e = NaN(1,2,3);
                        if norm(P1-P2) < Eps
                            tri_e(1,1,:) = element(e,[1,4,3]);
                        elseif norm(P1-P3) < Eps
                            tri_e(1,1,:) = element(e,[1,2,4]);
                        elseif norm(P2-P4) < Eps || norm(P3-P4) < Eps
                            tri_e(1,1,:) = element(e,[1,2,3]);
                        else
                            if norm(P2-P3) > norm(P1-P4)
                                tri_e(1,1,:) = element(e,[1,2,4]);
                                tri_e(1,2,:) = element(e,[1,4,3]);
                            else
                                tri_e(1,1,:) = element(e,[1,2,3]);
                                tri_e(1,2,:) = element(e,[2,4,3]);
                            end                                
                        end
                        tri(e,:,:) = tri_e;
                    end
                    tri = reshape(tri,size(tri,1)*size(tri,2),3);
                    tri(any(isnan(tri),2),:) = [];

                    %% Find h_max and store results
                    varCol{1}.h_max = max([norm2(P(tri(:,1),:)-P(tri(:,2),:)); 
                                        norm2(P(tri(:,1),:)-P(tri(:,3),:)); 
                                        norm2(P(tri(:,2),:)-P(tri(:,3),:))]);
                    varCol{1}.dofs = size(unique(tri,'rows','stable'),1);
                    varCol{1}.nepw = lambda(1)./varCol{1}.h_max;
                    varCol{1}.noElems = size(tri,1);
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     trisurf(tri,P(:,1),P(:,2),P(:,3), 'FaceColor', getColor(1))
%                     view(106,26) % sphere and cube
%                     axis off
%                     axis equal
%                     camlight
%                     ax = gca;               % get the current axis
%                     ax.Clipping = 'off';    % turn clipping off
%                     figureFullScreen(gcf)
% %                     
%                     export_fig(['../../graphics/sphericalShell/trianglesParm2_' num2str(varCol{1}.noElems)], '-png', '-transparent', '-r300')
%                     
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    
                    p_h = kirchApprTri(tri,P,v,varCol{1});
                case 'IGA'
                    varCol{1}.h_max = findMaxElementDiameter(varCol{1}.patches);
                    varCol{1}.nepw = varCol{1}.lambda./varCol{1}.h_max;
                    varCol{1}.dofs = varCol{1}.noDofs;
%                     p = calculateScatteredPressureBA(varCol{1}, Uc{1}, v, 0, plotFarField);
                    p_h = calculateScatteredPressureKDT(varCol{1}, v, plotFarField);
            end
        case 'RT'
            switch scatteringCase
                case 'MS'
                    d_vec = varCol{1}.d_vec;
                    p_h = zeros(size(d_vec,2),1);
%                     for i = 1:size(d_vec,2) %874%
                    plotFarField = task.plotFarField;
                    noIncDir = size(d_vec,2);
                    progressBars = varCol{1}.progressBars;
                    nProgressStepSize = ceil(noIncDir/1000);
                    if progressBars
                        ppm = ParforProgMon('Tracing rays: ', noIncDir, nProgressStepSize);
                    else
                        ppm = NaN;
                    end
%                     for i = 1:noIncDir
                    parfor i = 1:noIncDir
                        if progressBars && mod(i,nProgressStepSize) == 0
                            ppm.increment();
                        end
                        varColTemp2 = varCol{1};
                        varColTemp2.d_vec = d_vec(:,i);
%                         tic
                        varColTemp2 = createRays(varColTemp2);
%                         fprintf('\nCreating rays in %12f seconds.', toc)
%                         tic
                        varColTemp2 = traceRays(varColTemp2);    
%                         fprintf('\nTracing rays in %12f seconds.', toc)
%                         tic        
                        p_h(i) = calculateScatteredPressureRT(varColTemp2, v(i,:), plotFarField);
%                         fprintf('\nFar field in %12f seconds.', toc)
                    end
                otherwise
                    varCol{1} = createRays(varCol{1});
                    varCol{1} = traceRays(varCol{1});            
                    p_h = calculateScatteredPressureRT(varCol{1}, v, plotFarField);
            end
    end
    if ~runTasksInParallel
        fprintf('using %12f seconds.', toc)
    end
    task.results.p = p_h;
    task.results.abs_p = abs(p_h);
    task.results.TS = 20*log10(abs(p_h/P_inc));
    if analyticSolutionExist
        if plotFarField
            p_ref = varCol{1}.p_0(v);
%             p_ref = exactKDT(varCol{1}.k,varCol{1}.P_inc,parms.R_o);
        else
            p_ref = varCol{1}.p(v);
        end
        task.results.p_ref = p_ref;
        task.results.abs_p_ref = abs(p_ref);
        task.results.TS_ref = 20*log10(abs(p_ref/P_inc));

        task.results.error_pAbs = 100*abs(abs(p_ref)-abs(p_h))./abs(p_ref);
        task.results.error_p = 100*abs(p_ref-p_h)./abs(p_ref);
    end
end

% task.results.Uc{1} = Uc{1};
% if noDomains > 1
%     task.results.varCol{2} = varCol{2};
% %     task.results.Uc{2} = Uc{2};
% end
% if noDomains > 2
%     task.results.varCol{3} = varCol{3};
% %     task.results.Uc{3} = Uc{3};
% end
if strcmp(scatteringCase,'Ray')
    plotSolutionAlongRay
end

%% Calculate errors (if analyticSolutionExist) and plot result in Paraview
if ~task.useROM && ~strcmp(task.method,'RT')
    switch scatteringCase
        case {'BI', 'Sweep','Ray'}    
            tic
            if plotResidualError && analyticSolutionExist && strcmp(formulation,'GCBIE')
%                 close all
                fprintf(['\n%-' num2str(stringShift) 's'], 'Plotting zeros of residual ... ')
                plotGalerkinResidual(varCol{1},U);
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
                para.extraXiPts = eval(para.extraXiPts);
                para.extraEtaPts = eval(para.extraEtaPts);
                para.extraZetaPts = eval(para.extraZetaPts);
                testFun = @(v) -analytic(v) - P_inc(v);
                varCol{1}.testFun = testFun;
                if strcmp(task.method, 'KDT')
                    Uc{1} = zeros(varCol{1}.noDofs,1);
                end
        
                createParaviewFiles(varCol, 'U', Uc, 'para_options', para, 'e3Dss_options', e3Dss_options);

                if para.plotMesh
                    createVTKmeshFiles(varCol, 'U', Uc, 'para_options', para)
                end
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

if storeFullvarCol
    varColTemp = varCol{1};
else
    varColTemp.alpha = varCol{1}.alpha;
    varColTemp.beta = varCol{1}.beta;
    varColTemp.k = varCol{1}.k;
    varColTemp.f = varCol{1}.f;
    varColTemp.c_f = varCol{1}.c_f;
    if ~strcmp(task.method,'RT') && ~strcmp(task.method,'KDT')
        varColTemp.surfDofs = varCol{1}.surfDofs;
        varColTemp.dofs = varCol{1}.dofs;
        varColTemp.dofsAlg = (varCol{1}.dofs)^(1/3);
        varColTemp.h_max = varCol{1}.h_max;
        if isfield(varCol{1},'tau')
            varColTemp.tau = varCol{1}.tau;
        end
        varColTemp.nepw = varCol{1}.nepw;
        varColTemp.tot_time = varCol{1}.tot_time;
        varColTemp.noElems = varCol{1}.noElems;
        varColTemp.totNoElems = varCol{1}.totNoElems;
        if storeSolution
            varColTemp.U = varCol{1}.Uc;
        end
        if task.useROM
            varColTemp.U_sweep = U_sweep;
            varColTemp.U_sweep2 = U_sweep2;
        end
        if isfield(varCol{1},'N')
            varColTemp.N = varCol{1}.N;
        end
        if isfield(varCol{1},'timeBuildSystem')
            varColTemp.timeBuildSystem = varCol{1}.timeBuildSystem;
        end
        if isfield(varCol{1},'timeSolveSystem')
            varColTemp.timeSolveSystem = varCol{1}.timeSolveSystem;
        end
        if isfield(varCol{1},'totNoQP')
            varColTemp.totNoQP = varCol{1}.totNoQP;
        end
        if isfield(varCol{1},'totNoQPnonPolar')
            varColTemp.totNoQPnonPolar = varCol{1}.totNoQPnonPolar;
        end
    end 
    varColTemp.analyticSolutionExist = varCol{1}.analyticSolutionExist;
    varColTemp.boundaryMethod = varCol{1}.boundaryMethod;
end
task.varCol{1} = varColTemp;
