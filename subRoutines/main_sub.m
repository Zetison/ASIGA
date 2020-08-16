function task = main_sub(task,loopParameters,runTasksInParallel,resultsFolder)

stringShift = 40;
extractTaskData

defineFolderAndFilenames
if ~runTasksInParallel
    fprintf('\n\nRunning the case: %s', saveName)
end


setAndDefineParameters

if ~runTasksInParallel
    fprintf(['\n%-' num2str(stringShift) 's'], 'Extracting CAD data ... ')
    tic
end
varCol = createNURBSmesh(varCol, model, M, degree);
if ~runTasksInParallel
    fprintf('using %12f seconds.', toc)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Collect varables into varCol
collectVariables

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot geometry
% keyboard
if prePlot.plot3Dgeometry || prePlot.plot2Dgeometry
    tic
    fprintf(['\n%-' num2str(stringShift) 's'], 'Plotting geometry ... ')
    plotMeshAndGeometry
    fprintf('using %12f seconds.', toc)
end

%% Build connectivity
tic
if ~runTasksInParallel
    fprintf(['\n%-' num2str(stringShift) 's'], 'Generating IGA mesh ... ')
end
for i = 1:numel(varCol)
    varCol{i} = generateIGAmesh(varCol{i});
end
% if ~boundaryMethod && M(i_M) < 2
%     ratio = calcSA_Vratio(varCol{1});
%     SAVindex = L_gamma_a/2*ratio
% end

%% Find nodes to remove due to gluing and homogeneous dirichlet conditions
totNoElems = 0;
for i = 1:numel(varCol)
    varCol{i} = findDofsToRemove(varCol{i});
    totNoElems = totNoElems + varCol{i}.noElems;
end

% plotControlPts(fluid);
% dofsToRemove = varCol{1}.dofsToRemove;
% controlPts = varCol{1}.controlPts;
% plot3(controlPts(dofsToRemove,1),controlPts(dofsToRemove,2),controlPts(dofsToRemove,3),'o','color','blue','MarkerFaceColor','blue')

if ~runTasksInParallel
    fprintf('using %12f seconds.', toc)
    fprintf('\nTotal number of elements = %d', totNoElems)
end
varCol{1}.totNoElems = totNoElems;
if (prePlot.plot3Dgeometry || prePlot.plot2Dgeometry) && prePlot.abortAfterPlotting
    return
end
%% Build stiffness matrix
t_start = tic;
if strcmp(scatteringCase,'Sweep')
    U_sweep = cell(1,numel(f));
end

if strcmp(method,'RT') || strcmp(method,'KDT')
    omega = 2*pi*f;
    varCol{1}.omega = omega;
    varCol{1}.k = omega/varCol{1}.c_f;
    varCol{1}.lambda = 2*pi./varCol{1}.k;
    varCol{1}.f = f;
    getAnalyticSolutions
else
    for i_f = 1:numel(f)
        f_i = f(i_f);
        omega = 2*pi*f_i;
        for m = 1:numel(varCol)
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
        getAnalyticSolutions
        k = varCol{1}.k;
        rho = varCol{1}.rho;
        switch method
            case {'IE','ABC'}
                tic  
                if ~runTasksInParallel
                    fprintf(['\n%-' num2str(stringShift) 's'], 'Building outer fluid matrix ... ')
                end
                options = {'operator','Laplace',...
                           'fieldDimension', 1,...
                           'buildMassMatrix', 1};

                if strcmp(varCol{1}.coreMethod,'SEM')
                    [A_fluid_o, FF, varCol{1}.dofsToRemove] = buildSEMMatrices(varCol{1});
                    FF = 1/(rho(1)*omega^2)*FF;        
                    noDofs_tot = size(A_fluid_o,1);
                else
        %             [A_K, A_M] = buildGlobalMatrices(varCol{1}, options);
                    [A_K, A_M] = buildMatrices(varCol{1}, options);
                    A_fluid_o = A_K - k^2*A_M;
                end
                if clearGlobalMatrices && ~useROM
                    clear A_K A_M
                end
%                 nnzA_fluid_o = nnz(A_fluid_o);
                varCol{1}.timeBuildSystem = toc;
                if ~runTasksInParallel
                    fprintf('using %12f seconds.', toc)
                end

                if strcmp(method,'IE') && ~strcmp(varCol{1}.coreMethod,'SEM')
                    % Add contribution from infinite elements
                    tic
                    if ~runTasksInParallel
                        fprintf(['\n%-' num2str(stringShift) 's'], 'Building infinite element matrix ... ')
                    end
                    if useROM
%                         [A_gamma_a, A2_gamma_a, A3_gamma_a, newDofsToRemove] = addInfElements3_ROM(varCol{1});
                        [A_gamma_a, A2_gamma_a, A3_gamma_a, newDofsToRemove] = addInfElements3_ROM2(varCol{1});
                    else
                        [A_gamma_a, newDofsToRemove] = addInfElements3(varCol{1});
                    end

            %         [A_inf, newDofsToRemove] = addInfElements4(varCol{1}, k(1), Upsilon); 
                    varCol{1}.timeBuildSystem = varCol{1}.timeBuildSystem + toc;
                    if ~runTasksInParallel
                        fprintf('using %12f seconds.', toc)
                    end

                    noDofs_new = max(size(A_gamma_a));
                    varCol{1}.noDofs_new = noDofs_new;

                    if noDofs_new > varCol{1}.noDofs
                        A_fluid_o(noDofs_new, noDofs_new) = 0;
                        if useROM
                            A_K(noDofs_new, noDofs_new) = 0;
                            A_M(noDofs_new, noDofs_new) = 0;
                        end
                    end
                    dofsToRemove_old = varCol{1}.dofsToRemove;
                    varCol{1}.dofsToRemove = sort(unique([varCol{1}.dofsToRemove newDofsToRemove]));
                    varCol{1}.dofsToRemove_old = dofsToRemove_old;

                    noDofs_tot = noDofs_new;
                elseif strcmp(method,'ABC')
                    % Add contribution from infinite elements
                    tic
                    if ~runTasksInParallel
                        fprintf(['\n%-' num2str(stringShift) 's'], 'Building ABC matrix ... ')
                    end
                    A_gamma_a = addABC(varCol{1}); 

            %         [A_inf, newDofsToRemove] = addInfElements4(varCol{1}, k(1), Upsilon); 
                    varCol{1}.timeBuildSystem = varCol{1}.timeBuildSystem + toc;
                    if ~runTasksInParallel
                        fprintf('using %12f seconds.', toc)
                    end


                    noDofs_tot = varCol{1}.noDofs;
                end
                dofsToRemove = varCol{1}.dofsToRemove;
                if strcmp(varCol{1}.coreMethod,'SEM')
                    A_fluid_o = 1/(rho(1)*omega^2)*A_fluid_o; 
                else
                    if useROM
                        A_K = 1/(rho(1)*omega^2)*A_K; 
                        A_M = 1/(rho(1)*omega^2)*A_M; 
                        A_gamma_a = 1/(rho(1)*omega^2)*A_gamma_a; 
                        A2_gamma_a = 1/(rho(1)*omega^2)*A2_gamma_a; 
                        A3_gamma_a = 1/(rho(1)*omega^2)*A3_gamma_a; 
                    else
                        A_fluid_o = 1/(rho(1)*omega^2)*(A_fluid_o + A_gamma_a); 
                    end
                end

                if numel(varCol) > 1 
                    tic
                    % Solid matrices
                    if ~runTasksInParallel
                        fprintf(['\n%-' num2str(stringShift) 's'], 'Building solid matrix ... ')
                    end
                    options = {'operator','linearElasticity',...
                               'fieldDimension', 3,...
                               'buildMassMatrix', 1};
                    [A_K, A_M] = buildMatrices(varCol{2}, options);

                    A_solid = A_K-varCol{2}.rho*omega^2*A_M;
                    if clearGlobalMatrices
                        clear A_K A_M
                    end

                    noDofs_tot = noDofs_tot + varCol{2}.noDofs;
                    dofsToRemove = [varCol{2}.dofsToRemove (varCol{1}.dofsToRemove+varCol{2}.noDofs)];
                    varCol{1}.timeBuildSystem = varCol{1}.timeBuildSystem + toc;
                    if ~runTasksInParallel
                        fprintf('using %12f seconds.', toc)
                    end
                end

                if numel(varCol) > 2
                    tic        
                    % Inner fluid   
                    if ~runTasksInParallel
                        fprintf(['\n%-' num2str(stringShift) 's'], 'Building inner fluid matrix ... ')
                    end
                    options = {'operator','Laplace',...
                               'fieldDimension', 1,...
                               'buildMassMatrix', 1};
                    [A_K, A_M] = buildMatrices(varCol{3}, options);

                    A_fluid_i = 1/(varCol{3}.rho*omega^2)*(A_K - varCol{3}.k^2*A_M);
                    if clearGlobalMatrices
                        clear A_K A_M
                    end
                    if ~runTasksInParallel
                        fprintf('using %12f seconds.', toc)
                    end
                    noDofs_tot = noDofs_tot + varCol{3}.noDofs;
                    dofsToRemove = [varCol{3}.dofsToRemove (varCol{2}.dofsToRemove+varCol{3}.noDofs) (varCol{1}.dofsToRemove+varCol{3}.noDofs+varCol{2}.noDofs)];
                    varCol{1}.timeBuildSystem = varCol{1}.timeBuildSystem + toc;
                end


                % Collect all matrices
                A = sparse(noDofs_tot,noDofs_tot);
                switch numel(varCol)
                    case 1
                        varCol{1}.noDofs_tot = noDofs_tot;
                        A = A_fluid_o;   
                    case 2
                        varCol{2}.noDofs_tot = noDofs_tot;
                        A(1:varCol{2}.noDofs,1:varCol{2}.noDofs) = A_solid;  
                        A((varCol{2}.noDofs+1):noDofs_tot,(varCol{2}.noDofs+1):noDofs_tot) = A_fluid_o; 
                        shift = 0;
                    case 3
                        varCol{2}.noDofs_tot = noDofs_tot;
                        A(1:varCol{3}.noDofs,1:varCol{3}.noDofs) = A_fluid_i; 
                        A((varCol{3}.noDofs+1):(varCol{3}.noDofs+varCol{2}.noDofs),(varCol{3}.noDofs+1):(varCol{3}.noDofs+varCol{2}.noDofs)) = A_solid;
                        A((varCol{3}.noDofs+varCol{2}.noDofs+1):noDofs_tot,(varCol{3}.noDofs+varCol{2}.noDofs+1):noDofs_tot) = A_fluid_o;
                        shift = varCol{3}.noDofs;
                end

                % Apply coupling conditions    
                if numel(varCol) > 1 
                    tic
                    if ~runTasksInParallel
                        fprintf(['\n%-' num2str(stringShift) 's'], 'Building coupling matrix ... ')
                    end
%                     A_coupling = applyCouplingCondition_new(varCol{2},shift);
                    A_coupling = applyCouplingConditionPatches(varCol{2},varCol{1},varCol{2}.noDofs+shift,shift,noDofs_tot);
%                     A_coupling = applyCouplingConditionPatchesTest(varCol{2},varCol{1},shift,varCol{2}.noDofs,noDofs_tot);
%                     A_coupling = applyCouplingConditionPatches(varCol{2},varCol{1},shift,noDofs_tot);
                    A = A + A_coupling;
                    if ~runTasksInParallel
                        fprintf('using %12f seconds.', toc)
                    end
                    varCol{1}.timeBuildSystem = varCol{1}.timeBuildSystem + toc;
                end

                if numel(varCol) > 2 
                    tic
                    if ~runTasksInParallel
                        fprintf(['\n%-' num2str(stringShift) 's'], 'Building second coupling matrix ... ')
                    end
%                     A_coupling = applyCouplingCondition2(varCol{2},varCol{3},varCol{3}.noDofs);
                    A_coupling = -applyCouplingConditionPatches(varCol{3},varCol{2},shift,0,noDofs_tot);
                    A = A + A_coupling;
                    if ~runTasksInParallel
                        fprintf('using %12f seconds.', toc)
                    end
                    varCol{1}.timeBuildSystem = varCol{1}.timeBuildSystem + toc;
                end    

                if clearGlobalMatrices
                    clear A_fluid_o A_fluid_i A_solid A_coupling
                    if ~useROM
                        clear A_gamma_a
                    end
                end

                % Apply Neumann conditions
                tic
                if ~runTasksInParallel
                    fprintf(['\n%-' num2str(stringShift) 's'], 'Building right hand side vector ... ')
                end
                if numel(varCol) == 1
                    if ~strcmp(varCol{1}.coreMethod,'SEM')
                        if useROM
                            FF = 1/(rho(1)*omega^2)*applyHWBC_ROM(varCol{1},noVecs);
                        else
                            FF = 1/(rho(1)*omega^2)*applyHWBC(varCol{1},length(alpha_s));  
                        end
                    end
                else
                    varCol{2}.p_inc = p_inc;
                    varCol{2}.dp_inc = dp_inc;
%                     FF = applyNeumannCondition_CoupledProblem(varCol{2},omega,rho_f(1),length(alpha_s), shift);
                    FF = applyNeumannCondition_CoupledProblemPatches(varCol{2},varCol{1},omega,rho(1),length(alpha_s), shift);
                end 
                if ~runTasksInParallel
                    fprintf('using %12f seconds.', toc)
                end
                varCol{1}.timeBuildSystem = varCol{1}.timeBuildSystem + toc;

                %% Modify system of equations
                % Remove the rows and columns of the global matrix corresponding to
                % removed degrees of freedom
                if numel(varCol) > 1 
                    P1 = speye(size(A));
                    if numel(varCol) > 2
                        P1(1:varCol{3}.noDofs, 1:varCol{3}.noDofs) = ...
                            sqrt(omega^2*varCol{3}.rho)*speye(varCol{3}.noDofs);
                    end

                    P1(shift+1:shift+varCol{2}.noDofs, shift+1:shift+varCol{2}.noDofs) = ...
                        1/sqrt(max([omega^2*varCol{2}.rho, varCol{2}.C(1,1)]))*speye(varCol{2}.noDofs);
                    P1(shift+1+varCol{2}.noDofs:end, shift+1+varCol{2}.noDofs:end) = ...
                        sqrt(omega^2*rho(1))*speye(varCol{1}.noDofs_new);

                    A = P1*A*P1;
                    FF = P1*FF;

                    P1(dofsToRemove,:) = [];  
                    P1(:,dofsToRemove) = [];
                end

                if useROM
                    A_K(dofsToRemove,:) = [];  
                    A_K(:,dofsToRemove) = [];
                    A_M(dofsToRemove,:) = [];  
                    A_M(:,dofsToRemove) = [];
                    A_gamma_a(dofsToRemove,:) = [];  
                    A_gamma_a(:,dofsToRemove) = [];
                    A2_gamma_a(dofsToRemove,:) = [];  
                    A2_gamma_a(:,dofsToRemove) = [];
                    A3_gamma_a(dofsToRemove,:) = [];  
                    A3_gamma_a(:,dofsToRemove) = [];
                else
                    A(dofsToRemove,:) = [];  
                    A(:,dofsToRemove) = [];
                end

                FF(dofsToRemove,:) = [];
            case 'IENSG'
                tic
                if ~runTasksInParallel
                    fprintf(['\n%-' num2str(stringShift) 's'], 'Building infinite element matrix ... ')
                end
                [A, dofsToRemove] = infElementsNonSepGeom(varCol{1});  
                if ~runTasksInParallel
                    fprintf('using %12f seconds.', toc)
                end

                varCol{1}.dofsToRemove = dofsToRemove;
                noDofs_tot = varCol{1}.noDofs*varCol{1}.N;
                varCol{1}.noDofs_tot = noDofs_tot;

                % Apply Neumann conditions
                tic
                if ~runTasksInParallel
                    fprintf(['\n%-' num2str(stringShift) 's'], 'Building right hand side vector ... ')
                end

                FF = applyHWBC_nonSepGeom(varCol{1},alpha_s);

                varCol{1}.timeBuildSystem = toc;
                if ~runTasksInParallel
                    fprintf('using %12f seconds.', toc)
                end

                A(dofsToRemove,:) = [];  
                A(:,dofsToRemove) = [];

                FF(dofsToRemove,:) = [];
            case 'BEM'
                tic
                if ~runTasksInParallel
                    fprintf(['\n%-' num2str(stringShift) 's'], 'Building BEM matrix ... ')
                end
                dofsToRemove = varCol{1}.dofsToRemove;  
                switch formulation(1)
                    case 'G' % Galerkin
                        [A_fluid_o, FF_fluid_o, varCol{1}, C_mat] = buildGBEMmatrix(varCol{1},numel(varCol) > 1);  
                    case 'C' % Collocation
                        [A_fluid_o, FF_fluid_o, varCol{1}] = buildCBEMmatrix(varCol{1});  
                end
                noDofs_tot = varCol{1}.noDofs;
                if ~runTasksInParallel
                    fprintf('using %12f seconds.', toc)
                end
                if strcmp(formulation(end),'C')
                    tic
                    fprintf(['\n%-' num2str(stringShift) 's'], 'Building CHIEF matrix ... ')
                    [A_CHIEF, FF_CHIEF, varCol{1}] = buildCHIEFmatrix(varCol{1});
                    A_fluid_o = [A_fluid_o; A_CHIEF];
                    FF_fluid_o = [FF_fluid_o; FF_CHIEF];
                    if ~runTasksInParallel
                        fprintf('using %12f seconds.', toc)
                    end
                end
                varCol{1}.timeBuildSystem = toc(t_start);
                if numel(varCol) > 1 
                    tic
                    % Solid matrices
                    if ~runTasksInParallel
                        fprintf(['\n%-' num2str(stringShift) 's'], 'Building solid matrix ... ')
                    end
                    options = {'operator','linearElasticity',...
                               'fieldDimension', 3,...
                               'buildMassMatrix', 1};
                    [A_K, A_M] = buildMatrices(varCol{2}, options);

                    A_solid = A_K-varCol{2}.rho*omega^2*A_M;
                    if clearGlobalMatrices
                        clear A_K A_M
                    end

                    noDofs_tot = noDofs_tot + varCol{2}.noDofs;
                    dofsToRemove = [varCol{2}.dofsToRemove (varCol{1}.dofsToRemove+varCol{2}.noDofs)];
                    varCol{1}.timeBuildSystem = varCol{1}.timeBuildSystem + toc;
                    if ~runTasksInParallel
                        fprintf('using %12f seconds.', toc)
                    end
                end
                if numel(varCol) > 2
                    if ~runTasksInParallel
                        fprintf(['\n%-' num2str(stringShift) 's'], 'Building inner fluid matrix ... ')
                    end
                    dofsToRemove = [varCol{3}.dofsToRemove (varCol{2}.dofsToRemove+varCol{3}.noDofs) (varCol{1}.dofsToRemove+varCol{3}.noDofs+varCol{2}.noDofs)];

                    switch formulation(1)
                        case 'G' % Galerkin
                            [A_fluid_i, FF_fluid_i, varCol{3}, C_mat2_outer] = buildGBEMmatrix(varCol{3},numel(varCol) > 1);  
                        case 'C' % Collocation
                            error('Not implemented')
                    end
                    noDofs_tot = noDofs_tot + varCol{3}.noDofs;
                    if ~runTasksInParallel
                        fprintf('using %12f seconds.', toc)
                    end
                end
                % Collect all matrices
                A = sparse(noDofs_tot,noDofs_tot);
                FF = zeros(noDofs_tot,numel(varCol{1}.alpha_s));
                switch numel(varCol)
                    case 1
                        varCol{1}.noDofs_tot = noDofs_tot;
                        A = A_fluid_o;   
                        FF = FF_fluid_o;
                    case 2
                        varCol{1}.noDofs_tot = noDofs_tot;
                        varCol{2}.noDofs_tot = noDofs_tot;
                        A(1:varCol{2}.noDofs,1:varCol{2}.noDofs) = A_solid;  
                        A(varCol{2}.noDofs+1:end,(1:3*varCol{1}.noDofs)+varCol{2}.noDofs-3*varCol{1}.noDofs) ...
                                    = -rho(1)*omega^2*C_mat;  
                        A((varCol{2}.noDofs+1):noDofs_tot,(varCol{2}.noDofs+1):noDofs_tot) = A_fluid_o; 
                        shift = 0;
                        FF((varCol{2}.noDofs+1):noDofs_tot,:) = FF_fluid_o;
                    case 3
                        shift = varCol{3}.noDofs;
                        noDofsInner = varCol{3}.noDofsInner;
                        varCol{1}.noDofs_tot = noDofs_tot;
                        varCol{2}.noDofs_tot = noDofs_tot;
                        varCol{3}.noDofs_tot = noDofs_tot;
                        A(1:varCol{3}.noDofs,1:varCol{3}.noDofs) = A_fluid_i; 
                        A(noDofsInner+1:shift,shift+1:shift+3*(shift-noDofsInner)) = -varCol{3}.rho*omega^2*C_mat2_outer; 
                        A((varCol{3}.noDofs+1):(varCol{3}.noDofs+varCol{2}.noDofs),...
                          (varCol{3}.noDofs+1):(varCol{3}.noDofs+varCol{2}.noDofs)) = A_solid;
                        A(varCol{2}.noDofs+1+shift:end,(1:3*varCol{1}.noDofs)+varCol{2}.noDofs-3*varCol{1}.noDofs+shift) ...
                            = -rho(1)*omega^2*C_mat; 
                        A((varCol{3}.noDofs+varCol{2}.noDofs+1):noDofs_tot,...
                          (varCol{3}.noDofs+varCol{2}.noDofs+1):noDofs_tot) = A_fluid_o;
                        FF((shift+varCol{2}.noDofs+1):noDofs_tot,:) = FF_fluid_o;
                end
                % Apply coupling conditions    
                if numel(varCol) > 1 
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

                if numel(varCol) > 2 
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
                
                %% Modify system of equations
                % Remove the rows and columns of the global matrix corresponding to
                % removed degrees of freedom
                if numel(varCol) > 1 
                    P1 = speye(size(A));
%                     if numel(varCol) > 2
%                         P1(1:varCol{3}.noDofs, 1:varCol{3}.noDofs) = ...
%                             sqrt(omega^2*rho_f(2))*speye(varCol{3}.noDofs);
%                     end

                    P1(shift+1:shift+varCol{2}.noDofs, ...
                       shift+1:shift+varCol{2}.noDofs) = 1/max([omega^2*varCol{2}.rho, varCol{2}.C(1,1)])*speye(varCol{2}.noDofs);

%                     A = P1*A;
%                     FF = P1*FF;
                    A = A*P1;

                    P1(dofsToRemove,:) = [];  
                    P1(:,dofsToRemove) = [];
                end
%                 close all
%                 spy2(A)
                A(:,dofsToRemove) = [];
                if strcmp(formulation(1),'G')
                    A(dofsToRemove,:) = [];
                    FF(dofsToRemove,:) = [];
                end
%                 figure(2)
%                 spy(A)
%                 keyboard
            case 'BA'
                tic
                if ~runTasksInParallel
                    fprintf(['\n%-' num2str(stringShift) 's'], 'Building outer fluid matrix ... ')
                end
                [A_fluid_o, FF_fluid_o, varCol{1}] = buildBAmatrix(varCol{1});
                noDofs_tot = varCol{1}.noDofs;
                dofsToRemove = varCol{1}.dofsToRemove;
                varCol{1}.timeBuildSystem = toc;
                if ~runTasksInParallel
                    fprintf('using %12f seconds.', varCol{1}.timeBuildSystem)
                end
                if numel(varCol) > 1 
                    error('BA is not implemented in such a way that the displacement and pressure conditions at the interfaces is satisfied')
                    tic
                    % Solid matrix
                    if ~runTasksInParallel
                        fprintf(['\n%-' num2str(stringShift) 's'], 'Building solid matrix ... ')
                    end
                    varCol{2}.formulation = 'VL2E';
                    [A_solid, FF_solid, varCol{2}] = buildBAmatrix(varCol{2});
                    A_solid = kron(A_solid,eye(3));
                    FF_solid = FF_solid.';
                    FF_solid = FF_solid(:);
                    noDofs_tot = noDofs_tot + varCol{2}.noDofs;
                    dofsToRemove = [varCol{2}.dofsToRemove (varCol{1}.dofsToRemove+varCol{2}.noDofs)];
                    varCol{1}.timeBuildSystem = varCol{1}.timeBuildSystem + toc;
                    if ~runTasksInParallel
                        fprintf('using %12f seconds.', toc)
                    end
                end

                if numel(varCol) > 2
                    tic        
                    % Inner fluid   
                    if ~runTasksInParallel
                        fprintf(['\n%-' num2str(stringShift) 's'], 'Building inner fluid matrix ... ')
                    end
                    [A_fluid_i, FF_fluid_i, varCol{3}] = buildBAmatrix(varCol{3});

                    if ~runTasksInParallel
                        fprintf('using %12f seconds.', toc)
                    end
                    noDofs_tot = noDofs_tot + varCol{3}.noDofs;
                    dofsToRemove = [varCol{3}.dofsToRemove (varCol{2}.dofsToRemove+varCol{3}.noDofs) (varCol{1}.dofsToRemove+varCol{3}.noDofs+varCol{2}.noDofs)];
                    varCol{1}.timeBuildSystem = varCol{1}.timeBuildSystem + toc;
                end

                % Collect all matrices
                A = sparse(noDofs_tot,noDofs_tot);
                FF = zeros(noDofs_tot,1);
                switch numel(varCol)
                    case 1
                        varCol{1}.noDofs_tot = noDofs_tot;
                        A = A_fluid_o;   
                        FF = FF_fluid_o;
                    case 2
                        varCol{2}.noDofs_tot = noDofs_tot;
                        A(1:varCol{2}.noDofs,1:varCol{2}.noDofs) = A_solid;  
                        A((varCol{2}.noDofs+1):noDofs_tot,(varCol{2}.noDofs+1):noDofs_tot) = A_fluid_o; 
                        shift = 0;
                        FF(1:varCol{2}.noDofs) = FF_solid;
                        FF((varCol{2}.noDofs+1):noDofs_tot) = FF_fluid_o;
                    case 3
                        varCol{2}.noDofs_tot = noDofs_tot;
                        A(1:varCol{3}.noDofs,1:varCol{3}.noDofs) = A_fluid_i; 
                        A((varCol{3}.noDofs+1):(varCol{3}.noDofs+varCol{2}.noDofs),(varCol{3}.noDofs+1):(varCol{3}.noDofs+varCol{2}.noDofs)) = A_solid;
                        A((varCol{3}.noDofs+varCol{2}.noDofs+1):noDofs_tot,(varCol{3}.noDofs+varCol{2}.noDofs+1):noDofs_tot) = A_fluid_o;
                        shift = varCol{3}.noDofs;
                        FF(1:varCol{3}.noDofs) = FF_fluid_i;
                        FF((varCol{3}.noDofs+1):(varCol{3}.noDofs+varCol{2}.noDofs)) = FF_solid;
                        FF((varCol{3}.noDofs+varCol{2}.noDofs+1):noDofs_tot) = FF_fluid_o;
                end

                %% Modify system of equations
                % Remove the rows and columns of the global matrix corresponding to
                % removed degrees of freedom
                if numel(varCol) > 1 
                    P1 = speye(size(A));

                    P1(dofsToRemove,:) = [];  
                    P1(:,dofsToRemove) = [];
                end

                A(dofsToRemove,:) = [];  
                A(:,dofsToRemove) = [];

                FF(dofsToRemove,:) = [];

                varCol{1}.timeBuildSystem = toc;
            case 'MFS'
                tic
                if ~runTasksInParallel
                    fprintf(['\n%-' num2str(stringShift) 's'], 'Building MFS matrix ... ')
                end
                dofsToRemove = varCol{1}.dofsToRemove;  
                [A, FF, varCol{1}] = buildMFSmatrix(varCol{1});  
    %             rmsA = rms(rms(A))
    %             condA = cond(A)
                noDofs_tot = varCol{1}.noDofs;
                varCol{1}.timeBuildSystem = toc;
                if ~runTasksInParallel
                    fprintf('using %12f seconds.', varCol{1}.timeBuildSystem)
                end
            case 'KDT'
                tic
                if ~runTasksInParallel
                    fprintf(['\n%-' num2str(stringShift) 's'], 'Building outer fluid matrix ... ')
                end
                [A, FF, varCol{1}] = bestApproximationVec(varCol{1});
                dofsToRemove = varCol{1}.dofsToRemove;
                A(dofsToRemove,:) = [];  
                A(:,dofsToRemove) = [];

                FF(dofsToRemove,:) = [];
                noDofs_tot = varCol{1}.noDofs;
                varCol{1}.timeBuildSystem = toc;
                if ~runTasksInParallel
                    fprintf('using %12f seconds.', varCol{1}.timeBuildSystem)
                end
        end
        switch method
            case {'IE','IENSG','BEM','BA','ABC','MFS'}
                if ~runTasksInParallel
                    fprintf(['\n%-' num2str(stringShift) 's'], 'Total time building system ... ')
                    fprintf(' was  %12f seconds.', varCol{1}.timeBuildSystem)
                end
                % condA = condest(A)
                %% SOLVE SYSTEM
                if 0 %strcmp(BC,'SHBC') && (strcmp(method,'IE') || strcmp(method,'ABC')) && M > 6 && ~strcmp(coreMethod,'IGA') && ~(strcmp(coreMethod,'linear_FEM') && M == 7)
                %     solver = 'cg';
                    solver = 'gmres';
                    % solver = 'cgs';
                %     solver = 'LU';
                else
                    solver = 'LU';
                end
                if strcmp(BC,'SHBC')
                    if useROM && ~strcmp(method,'BA')
                        A_M = rho(1)*omega^2*A_M;
                        A_K = rho(1)*omega^2*A_K;
                        A_gamma_a = rho(1)*omega^2*A_gamma_a;
                        A2_gamma_a = rho(1)*omega^2*A2_gamma_a;
                        A3_gamma_a = rho(1)*omega^2*A3_gamma_a;
                        FF = rho(1)*omega^2*FF;
                    else
                        A = rho(1)*omega^2*A;
                        FF = rho(1)*omega^2*FF;
                    end
                end
                if ~strcmp(solver,'LU')
                    tic
                    preconditioner = 'ilu';
                %     preconditioner = 'SSOR';
                    if ~runTasksInParallel
                        fprintf(['\n%-' num2str(stringShift) 's'], ['Creating preconditioner (' preconditioner ') ... '])
                    end
                    switch preconditioner
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
                switch solver
                    case 'cg'
                %         [UU,~,~,it1,rv1] = bicgstab(A,FF,1e-20,1000,L_A,U_A);
                %         [UU,~,~,it1,rv1] = bicgstabl(A,FF,1e-20,1000,L_A,U_A);
                        [UU,~,~,it1,rv1] = lsqr(A,FF,1e-20,1000,L_A,U_A);
                %         [UU,~,~,it1,rv1] = bicg(A,FF,1e-20,1000,L_A,U_A);
                %         [UU,~,~,it1,rv1] = cgs(A,FF,1e-20,1000,L_A,U_A);
                        if clearGlobalMatrices
                            clear L_A U_A
                        end
                    case 'gmres'
                        noRestarts = []; % change this to save memory
                        % noRestarts = 2;
                        [UU,~,~,it1,rv1] = gmres(A,FF,noRestarts,1e-20,1000,L_A,U_A);
                        if clearGlobalMatrices
                            clear L_A U_A
                        end
                    case 'LU'
                        if useROM
                            switch method
                                case 'BA'
                                    UU = A\FF;
                                otherwise
%                                     A = A_K - k_1^2*A_M + A_gamma_a;
%                                     A2 = - 2*k_1*A_M + A2_gamma_a;
%                                     A3 = - 2*A_M + A3_gamma_a;
                                    A = A_K - k_1^2*A_M + k_1^2*A_gamma_a + k_1*A2_gamma_a + A3_gamma_a;
                                    A2 = - 2*k_1*A_M + 2*k_1*A_gamma_a + A2_gamma_a;
                                    A3 = - 2*A_M + 2*A_gamma_a;
                                    if task.useDGP
                                        varCol{1}.A_K = A_K;
                                        varCol{1}.A_M = A_M;
                                        varCol{1}.A_gamma_a = A_gamma_a;
                                        varCol{1}.A2_gamma_a = A2_gamma_a;
                                        varCol{1}.A3_gamma_a = A3_gamma_a;
                                    end
                                    UU = zeros(size(FF));
                %                     [L_A,U_A,P] = lu(A);
                                    dA = decomposition(A,'lu');
            %                         fprintf('using %12f seconds.', toc)
            %                         fprintf(['\n%-' num2str(stringShift) 's'], 'Computing ROM solution ... ')
                                    for i = 1:noVecs
                                        j = i-1;
                                        b = FF(:,i);
                                        if j > 0
                                            b = b - j*A2*UU(:,i-1);
                                        end
                                        if j > 1
                                            b = b - j*(j-1)/2*A3*UU(:,i-2);
                                        end
                %                         UU(:,i) = U_A\(L_A\(P*b));
                                        UU(:,i) = dA\b;
                %                         i
                                    end
                            end
    %                         fprintf('using %12f seconds.', toc)
                        else
                            if strcmp(method,'MFS') && strcmp(formulation,'SS')
                                Pinv = diag(1./max(abs(A)));
                                UU = diag(Pinv).*((A*Pinv)\FF);
                            else
                                UU = A\FF;
                            end
                        end
                end
                if ~runTasksInParallel
                    fprintf('using %12f seconds.', toc)
                end
                varCol{1}.timeSolveSystem = toc;
                % it1(end),rv1(end)
                % UU2 = UU;
                % load('UU');
                % norm(UU-UU2)/norm(UU)
                % keyboard
                % normUU = norm(UU)
                % condest(A)
                % keyboard
                if computeCondNumber && (size(A,1) == size(A,2))
                    rng('default') % for reproducibility in condest
                    condNumber = condest(A);
        %             condNumber = cond(full(A))
                else
                    condNumber = NaN;
                end
                % condNumber
                task.results.cond_number = condNumber;
                % condest(A)
                % [L_A,U_A] = ilu(A,struct('type','nofill'));
                % [L_A,U_A] = ilu(A,struct('type','crout','milu','row','droptol',0.1));
                % condest(A)
                % issymmetric2(A,1e-10)
                % keyboard
                % condSparse = condest(A)
                % condFull = cond(full(A))
                % keyboard

                % switch solver
                %     case 'cg'
                %         figure(47)
                %         semilogy(0:length(rv1)-1,rv1,'-o');
                %         xlabel('Iteration number');
                %         ylabel('Residual');
                %     case 'gmres'
                %         figure(47)
                %         semilogy(0:length(rv1)-1,rv1,'-o');
                %         xlabel('Iteration number');
                %         ylabel('Residual');
                % end
                % savefig('residualErrors')
                % figure(48)
                % spy2(A)
                % condest((L_A*U_A)\A)
                % tic; UU2 = A\FF; toc
                % max(abs(UU-UU2)./abs(UU2))
                % return

                if numel(varCol) > 1 && ~strcmp(method,'BA')
                    UU = P1*UU;
                end

                actualNoDofs = size(A,1);
                if ~runTasksInParallel
                    fprintf('\nNumber of degrees of freedom = %d', actualNoDofs)
                end
                % fprintf('\nMemory ratio = %f', ((fluid.degree(1)+1)^6*varCol{1}.noElems)/nnz(A_fluid_o))
                varCol{1}.actualNoDofs = actualNoDofs;
                % nXn = varCol{1}.nurbs.number(1);
                % nEn = varCol{1}.nurbs.number(2);
                % nZn = varCol{1}.nurbs.number(3);
                % nXn*nEn*(nZn+N-1) - (nZn+N-1)*(nEn-2) - 2*(nXn-1)*(nZn+N-1)
                if clearGlobalMatrices
                    clear A FF P1
                end
                postProcessSolution
                if boundaryMethod
                    varCol{1}.tau = computeTau(varCol{1});
                end
        end
        if useROM && strcmp(scatteringCase,'Sweep')
            U_sweep{i_f} = Uc{1}(1:varCol{1}.noDofs,:);
            U_sweep2{i_f} = Uc{1};
            if strcmp(BC,'SSBC') || strcmp(BC,'NNBC')
                error('not implemented due to noDofs')
            end
        end
        if ~useROM
            calculateErrors
        end
        if strcmp(scatteringCase,'Sweep') && size(k,2) > 1
            fprintf('\nTotal time spent on frequency: %12f\n', toc(t_freq))  
        end
    end
end
if useROM
    varCol{1}.U_sweep = U_sweep;
    varCol{1}.U_sweep2 = U_sweep2;
    varCol{1}.k = k;
end

%% Compute scattered pressure    
tic
if calculateFarFieldPattern && ~useROM
    if ~runTasksInParallel
        fprintf(['\n%-' num2str(stringShift) 's'], 'Computing far-field pattern ... ')
    end
    v = getFarFieldPoints(task.alpha,task.beta,task.r);

    switch method
        case {'IE','ABC','IENSG','BA','BEM'}
            p_h = calculateScatteredPressure(varCol, Uc, v, 0, plotFarField);
        case 'MFS'
            p_h = calculateScatteredPressureMFS(varCol{1}, Uc{1}, v, plotFarField);
        case 'KDT'
            switch coreMethod
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
% if numel(varCol) > 1
%     task.results.varCol{2} = varCol{2};
% %     task.results.Uc{2} = Uc{2};
% end
% if numel(varCol) > 2
%     task.results.varCol{3} = varCol{3};
% %     task.results.Uc{3} = Uc{3};
% end
if strcmp(scatteringCase,'Ray')
    plotSolutionAlongRay
end

%% Calculate errors (if analyticSolutionExist) and plot result in Paraview
if ~useROM && ~strcmp(method,'RT')
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
                savefig([resultsFolder '/' saveName])
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
    %             savefig([resultsFolderName '/' saveName '_surfPlot_mesh' num2str(M) '_formulation_' formulation '_degree' num2str(max(fluid.degree)) '.fig'])
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
                vtfFileName = [resultsFolderNameParaview '/' saveName];
                if isempty(para.name)
                    para.name = vtfFileName;
                end


                para.plotArtificialBoundary = para.plotArtificialBoundary && (strcmp(method,'IE') || strcmp(method,'ABC') || strcmp(method,'PML'));
                para.extraXiPts = eval(para.extraXiPts);
                para.extraEtaPts = eval(para.extraEtaPts);
                para.extraZetaPts = eval(para.extraZetaPts);
                testFun = @(v) -analytic(v) - P_inc(v);
                varCol{1}.testFun = testFun;
                if strcmp(method, 'KDT')
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

if storeFullvarCol || useROM
    varColTemp = varCol{1};
else
    varColTemp.alpha = varCol{1}.alpha;
    varColTemp.beta = varCol{1}.beta;
    varColTemp.k = varCol{1}.k;
    varColTemp.f = varCol{1}.f;
    varColTemp.c_f = varCol{1}.c_f;
    if ~strcmp(method,'RT') && ~strcmp(method,'KDT')
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
            varColTemp.U = varCol{1}.U;
        end
        if useROM
            varColTemp.U_sweep = U_sweep;
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
