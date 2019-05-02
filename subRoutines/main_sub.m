function task = main_sub(task,loopParameters,runTasksInParallel)

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
[varCol, fluid, solid, fluid_i] = createNURBSmesh(varCol, parms, model, M, degree);
if ~runTasksInParallel
    fprintf('using %12f seconds.', toc)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Collect varables into varCol, varCol_solid and varCol_fluid_i
collectVariables

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot geometry
% keyboard
if plot3Dgeometry || plot2Dgeometry
    tic
    fprintf(['\n%-' num2str(stringShift) 's'], 'Plotting geometry ... ')
    plotMeshAndGeometry
    fprintf('using %12f seconds.', toc)
%     return
end

%% Build connectivity
tic
if ~runTasksInParallel
    fprintf(['\n%-' num2str(stringShift) 's'], 'Generating IGA mesh ... ')
end
if varCol.boundaryMethod
    varCol = generateIGA2DMesh_new(varCol);
else
    varCol = generateIGA3DMesh_new(varCol);
end

if useSolidDomain   
    varCol_solid = generateIGA3DMesh_new(varCol_solid);
end

if useInnerFluidDomain  
    varCol_fluid_i = generateIGA3DMesh_new(varCol_fluid_i);
end
% if ~boundaryMethod && M(i_M) < 2
%     ratio = calcSA_Vratio(varCol);
%     SAVindex = L_gamma_a/2*ratio
% end

%% Find nodes to remove due to gluing and homogeneous dirichlet conditions
varCol = findDofsToRemove(varCol);

% plotControlPts(fluid);
% dofsToRemove = varCol.dofsToRemove;
% controlPts = varCol.controlPts;
% plot3(controlPts(dofsToRemove,1),controlPts(dofsToRemove,2),controlPts(dofsToRemove,3),'o','color','blue','MarkerFaceColor','blue')

totNoElems = varCol.noElems;
if useSolidDomain  
    varCol_solid = findDofsToRemove(varCol_solid);
    totNoElems = totNoElems + varCol_solid.noElems;
end  

if useInnerFluidDomain   
    varCol_fluid_i = findDofsToRemove(varCol_fluid_i);
    totNoElems = totNoElems + varCol_fluid_i.noElems;
end
if ~runTasksInParallel
    fprintf('using %12f seconds.', toc)
    fprintf('\nTotal number of elements = %d', totNoElems)
end
varCol.totNoElems = totNoElems;
if plot3Dgeometry || plot2Dgeometry
    return
end
%% Build stiffness matrix
t_start = tic;
if strcmp(scatteringCase,'Sweep')
    U_sweep = cell(1,numel(k));
end

for i_k = 1:size(k,2)
    t_freq = tic;
    k_1 = k(1,i_k);
    varCol.k = k_1;
    omega = k_1*varCol.c_f(1);
    varCol.omega = omega;
    varCol.f = omega/(2*pi);
    getAnalyticSolutions
    switch method
        case {'IE','ABC'}
            tic  
            if ~runTasksInParallel
                fprintf(['\n%-' num2str(stringShift) 's'], 'Building outer fluid matrix ... ')
            end
            options = {'operator','Laplace',...
                       'fieldDimension', 1,...
                       'buildMassMatrix', 1};

            if strcmp(varCol.coreMethod,'SEM')
                [A_fluid_o, FF, varCol.dofsToRemove] = buildSEMMatrices(varCol);
                FF = 1/(rho_f(1)*omega^2)*FF;        
                noDofs_tot = size(A_fluid_o,1);
            else
    %             [A_K, A_M] = buildGlobalMatrices(varCol, options);
                [A_K, A_M] = buildGlobalMatricesVec(varCol, options);
                A_fluid_o = A_K - k_1^2*A_M;
            end
            if clearGlobalMatrices && ~useROM
                clear A_K A_M
            end
            nnzA_fluid_o = nnz(A_fluid_o);
            varCol.timeBuildSystem = toc;
            if ~runTasksInParallel
                fprintf('using %12f seconds.', toc)
            end

            if strcmp(method,'IE') && ~strcmp(varCol.coreMethod,'SEM')
                % Add contribution from infinite elements
                tic
                if ~runTasksInParallel
                    fprintf(['\n%-' num2str(stringShift) 's'], 'Building infinite element matrix ... ')
                end
                if useROM
                    [A_gamma_a, A2_gamma_a, A3_gamma_a, newDofsToRemove] = addInfElements3_ROM(varCol);
                else
                    [A_gamma_a, newDofsToRemove] = addInfElements3(varCol);
                end

        %         [A_inf, newDofsToRemove] = addInfElements4(varCol, k(1), Upsilon); 
                varCol.timeBuildSystem = varCol.timeBuildSystem + toc;
                if ~runTasksInParallel
                    fprintf('using %12f seconds.', toc)
                end

                noDofs_new = max(size(A_gamma_a));
                varCol.noDofs_new = noDofs_new;

                if noDofs_new > varCol.noDofs
                    A_fluid_o(noDofs_new, noDofs_new) = 0;
                    if useROM
                        A_K(noDofs_new, noDofs_new) = 0;
                        A_M(noDofs_new, noDofs_new) = 0;
                    end
                end
                dofsToRemove_old = varCol.dofsToRemove;
                varCol.dofsToRemove = sort(unique([varCol.dofsToRemove newDofsToRemove]));
                varCol.dofsToRemove_old = dofsToRemove_old;

                noDofs_tot = noDofs_new;
            elseif strcmp(method,'ABC')
                % Add contribution from infinite elements
                tic
                if ~runTasksInParallel
                    fprintf(['\n%-' num2str(stringShift) 's'], 'Building ABC matrix ... ')
                end
                A_gamma_a = addABC(varCol); 

        %         [A_inf, newDofsToRemove] = addInfElements4(varCol, k(1), Upsilon); 
                varCol.timeBuildSystem = varCol.timeBuildSystem + toc;
                if ~runTasksInParallel
                    fprintf('using %12f seconds.', toc)
                end


                noDofs_tot = varCol.noDofs;
            end
            dofsToRemove = varCol.dofsToRemove;
            if strcmp(varCol.coreMethod,'SEM')
                A_fluid_o = 1/(rho_f(1)*omega^2)*A_fluid_o; 
            else
                if useROM
                    A_K = 1/(rho_f(1)*omega^2)*A_K; 
                    A_M = 1/(rho_f(1)*omega^2)*A_M; 
                    A_gamma_a = 1/(rho_f(1)*omega^2)*A_gamma_a; 
                    A2_gamma_a = 1/(rho_f(1)*omega^2)*A2_gamma_a; 
                    A3_gamma_a = 1/(rho_f(1)*omega^2)*A3_gamma_a; 
                else
                    A_fluid_o = 1/(rho_f(1)*omega^2)*(A_fluid_o + A_gamma_a); 
                end
            end

            if useSolidDomain 
                tic
                % Solid matrices
                if ~runTasksInParallel
                    fprintf(['\n%-' num2str(stringShift) 's'], 'Building solid matrix ... ')
                end
                options = {'operator','linearElasticity',...
                           'fieldDimension', 3,...
                           'buildMassMatrix', 1};
                [A_K, A_M] = buildGlobalMatrices(varCol_solid, options);

                A_solid = A_K-rho_s*omega^2*A_M;
                if clearGlobalMatrices
                    clear A_K A_M
                end

                noDofs_tot = noDofs_tot + varCol_solid.noDofs;
                dofsToRemove = [varCol_solid.dofsToRemove (varCol.dofsToRemove+varCol_solid.noDofs)];
                varCol.timeBuildSystem = varCol.timeBuildSystem + toc;
                if ~runTasksInParallel
                    fprintf('using %12f seconds.', toc)
                end
            end

            if useInnerFluidDomain
                tic        
                % Inner fluid   
                if ~runTasksInParallel
                    fprintf(['\n%-' num2str(stringShift) 's'], 'Building inner fluid matrix ... ')
                end
                options = {'operator','laplace',...
                           'fieldDimension', 1,...
                           'buildMassMatrix', 1};
                [A_K, A_M] = buildGlobalMatrices(varCol_fluid_i, options);

                k_2 = k(2,i_k);
                A_fluid_i = 1/(rho_f(2)*omega^2)*(A_K - k_2^2*A_M);
                if clearGlobalMatrices
                    clear A_K A_M
                end
                if ~runTasksInParallel
                    fprintf('using %12f seconds.', toc)
                end
                noDofs_tot = noDofs_tot + varCol_fluid_i.noDofs;
                dofsToRemove = [varCol_fluid_i.dofsToRemove (varCol_solid.dofsToRemove+varCol_fluid_i.noDofs) (varCol.dofsToRemove+varCol_fluid_i.noDofs+varCol_solid.noDofs)];
                varCol.timeBuildSystem = varCol.timeBuildSystem + toc;
            end


            % Collect all matrices
            A = sparse(noDofs_tot,noDofs_tot);
            if ~useSolidDomain
                varCol.noDofs_tot = noDofs_tot;
                A = A_fluid_o;   
            elseif useSolidDomain && ~useInnerFluidDomain
                varCol_solid.noDofs_tot = noDofs_tot;
                A(1:varCol_solid.noDofs,1:varCol_solid.noDofs) = A_solid;  
                A((varCol_solid.noDofs+1):noDofs_tot,(varCol_solid.noDofs+1):noDofs_tot) = A_fluid_o; 
                shift = 0;
            else 
                varCol_solid.noDofs_tot = noDofs_tot;
                A(1:varCol_fluid_i.noDofs,1:varCol_fluid_i.noDofs) = A_fluid_i; 
                A((varCol_fluid_i.noDofs+1):(varCol_fluid_i.noDofs+varCol_solid.noDofs),(varCol_fluid_i.noDofs+1):(varCol_fluid_i.noDofs+varCol_solid.noDofs)) = A_solid;
                A((varCol_fluid_i.noDofs+varCol_solid.noDofs+1):noDofs_tot,(varCol_fluid_i.noDofs+varCol_solid.noDofs+1):noDofs_tot) = A_fluid_o;
                shift = varCol_fluid_i.noDofs;
            end

            % Apply coupling conditions    
            if useSolidDomain 
                tic
                if ~runTasksInParallel
                    fprintf(['\n%-' num2str(stringShift) 's'], 'Building coupling matrix ... ')
                end
                A_coupling = applyCouplingCondition_new(varCol_solid,shift);
                A = A + A_coupling;
                if ~runTasksInParallel
                    fprintf('using %12f seconds.', toc)
                end
                varCol.timeBuildSystem = varCol.timeBuildSystem + toc;
            end

            if useInnerFluidDomain 
                tic
                if ~runTasksInParallel
                    fprintf(['\n%-' num2str(stringShift) 's'], 'Building second coupling matrix ... ')
                end
                A_coupling = applyCouplingCondition2(varCol_solid,varCol_fluid_i,varCol_fluid_i.noDofs);
                A = A + A_coupling;
                if ~runTasksInParallel
                    fprintf('using %12f seconds.', toc)
                end
                varCol.timeBuildSystem = varCol.timeBuildSystem + toc;
            end    

            if clearGlobalMatrices
                clear A_fluid_o A_fluid_i A_solid A_coupling
            end

            % Apply Neumann conditions
            tic
            if ~runTasksInParallel
                fprintf(['\n%-' num2str(stringShift) 's'], 'Building right hand side vector ... ')
            end
            if ~useSolidDomain
                if ~strcmp(varCol.coreMethod,'SEM')
                    if useROM
                        FF = 1/(rho_f(1)*omega^2)*applyHWBC_ROM(varCol,noVecs);
                    else
                        FF = 1/(rho_f(1)*omega^2)*applyHWBC(varCol,length(alpha_s));  
                    end
                end
            else
                varCol_solid.p_inc = p_inc;
                varCol_solid.dp_inc = dp_inc;
                FF = applyNeumannCondition_CoupledProblem(varCol_solid,omega,rho_f(1),length(alpha_s), shift);
            end 
            if ~runTasksInParallel
                fprintf('using %12f seconds.', toc)
            end
            varCol.timeBuildSystem = varCol.timeBuildSystem + toc;

            %% Modify system of equations
            % Remove the rows and columns of the global matrix corresponding to
            % removed degrees of freedom
            if useSolidDomain 
                P1 = speye(size(A));
                if useInnerFluidDomain
                    P1(1:varCol_fluid_i.noDofs, 1:varCol_fluid_i.noDofs) = ...
                        sqrt(omega^2*rho_f(2))*speye(varCol_fluid_i.noDofs);
                end

                P1(shift+1:shift+varCol_solid.noDofs, shift+1:shift+varCol_solid.noDofs) = ...
                    1/sqrt(max([omega^2*rho_s, C(1,1)]))*speye(varCol_solid.noDofs);
                P1(shift+1+varCol_solid.noDofs:end, shift+1+varCol_solid.noDofs:end) = ...
                    sqrt(omega^2*rho_f(1))*speye(varCol.noDofs_new);

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
            [A, dofsToRemove] = infElementsNonSepGeom(varCol);  
            if ~runTasksInParallel
                fprintf('using %12f seconds.', toc)
            end

            varCol.dofsToRemove = dofsToRemove;
            noDofs_tot = varCol.noDofs*varCol.N;
            varCol.noDofs_tot = noDofs_tot;

            % Apply Neumann conditions
            tic
            if ~runTasksInParallel
                fprintf(['\n%-' num2str(stringShift) 's'], 'Building right hand side vector ... ')
            end

            FF = applyHWBC_nonSepGeom(varCol,alpha_s);

            varCol.timeBuildSystem = toc;
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
            dofsToRemove = varCol.dofsToRemove;  
            switch formulation(1)
                case 'G' % Galerkin
                    switch formulation(2)
                        case 'R' % Regularized
%                             [A, FF] = buildRBEMmatrix_galerkin(varCol);  
                            [A, FF] = buildRBEMmatrix_galerkinVec(varCol);  
                        otherwise
%                             [A, FF] = buildBEMmatrix_galerkin(varCol);  
%                             [A, FF] = buildBEMmatrix_galerkinVec(varCol);  
                            [A, FF] = buildBEMmatrix_galerkinVec2(varCol);  
%                             [A, FF] = buildBEMmatrix_galerkinVec3(varCol);  
                    end
                    A(dofsToRemove,:) = [];
                    FF(dofsToRemove,:) = [];
                case 'C' % Collocation
                    switch formulation(2)
                        case 'R' % Regularized
                            [A, FF] = buildRBEMmatrixVec(varCol);  
%                             [A, FF] = buildRBEMmatrix(varCol);  
                        otherwise
%                             [A, FF] = buildBEMmatrixVec(varCol);  
                            [A, FF] = buildBEMmatrixVec2(varCol);  
%                             [A, FF] = buildBEMmatrixVec3(varCol);  
    %                         [A, FF] = buildBEMmatrix(varCol);  
                    end
            end 

            A(:,dofsToRemove) = [];
            noDofs_tot = varCol.noDofs;
            varCol.timeBuildSystem = toc;
            if ~runTasksInParallel
                fprintf('using %12f seconds.', varCol.timeBuildSystem)
            end
        case 'BA'
            tic
            if ~runTasksInParallel
                fprintf(['\n%-' num2str(stringShift) 's'], 'Building BA matrix ... ')
            end
            [A, FF] = bestApproximationVec(varCol);

            dofsToRemove = varCol.dofsToRemove;  
            noDofs_tot = varCol.noDofs;

            A(dofsToRemove,:) = [];  
            A(:,dofsToRemove) = [];

            FF(dofsToRemove,:) = [];

            varCol.timeBuildSystem = toc;
            if ~runTasksInParallel
                fprintf('using %12f seconds.', varCol.timeBuildSystem)
            end
        case 'MFS'
            tic
            if ~runTasksInParallel
                fprintf(['\n%-' num2str(stringShift) 's'], 'Building MFS matrix ... ')
            end
            dofsToRemove = varCol.dofsToRemove;  
            [A, FF, varCol] = buildMFSmatrix(varCol);  
            rmsA = rms(rms(A))
            condA = cond(A)
            noDofs_tot = varCol.noDofs;
            varCol.timeBuildSystem = toc;
            if ~runTasksInParallel
                fprintf('using %12f seconds.', varCol.timeBuildSystem)
            end
    end
    switch method
        case {'IE','IENSG','BEM','BA','ABC','MFS'}
            if ~runTasksInParallel
                fprintf(['\n%-' num2str(stringShift) 's'], 'Total time building system ... ')
                fprintf(' was  %12f seconds.', varCol.timeBuildSystem)
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
                    A_M = rho_f(1)*omega^2*A_M;
                    A_K = rho_f(1)*omega^2*A_K;
                    A_gamma_a = rho_f(1)*omega^2*A_gamma_a;
                    A2_gamma_a = rho_f(1)*omega^2*A2_gamma_a;
                    A3_gamma_a = rho_f(1)*omega^2*A3_gamma_a;
                    FF = rho_f(1)*omega^2*FF;
                else
                    A = rho_f(1)*omega^2*A;
                    FF = rho_f(1)*omega^2*FF;
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
                                A = A_K - k_1^2*A_M + A_gamma_a;
                                A2 = - 2*k_1*A_M + A2_gamma_a;
                                A3 = - 2*A_M + A3_gamma_a;
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
                        UU = A\FF;
                    end
            end
            if ~runTasksInParallel
                fprintf('using %12f seconds.', toc)
            end
            varCol.timeSolveSystem = toc;
            % it1(end),rv1(end)
            % UU2 = UU;
            % load('UU');
            % norm(UU-UU2)/norm(UU)
            % keyboard
            % normUU = norm(UU)
            % condest(A)
            % keyboard
            if computeCondNumber
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

            if useSolidDomain 
                UU = P1*UU;
            end

            actualNoDofs = size(A,1);
            if ~runTasksInParallel
                fprintf('\nNumber of degrees of freedom = %d', actualNoDofs)
            end
            % fprintf('\nMemory ratio = %f', ((fluid.degree(1)+1)^6*varCol.noElems)/nnz(A_fluid_o))
            varCol.actualNoDofs = actualNoDofs;
            % nXn = varCol.nurbs.number(1);
            % nEn = varCol.nurbs.number(2);
            % nZn = varCol.nurbs.number(3);
            % nXn*nEn*(nZn+N-1) - (nZn+N-1)*(nEn-2) - 2*(nXn-1)*(nZn+N-1)
            if clearGlobalMatrices
                clear A FF P1
            end
            postProcessSolution
    end
    if strcmp(scatteringCase,'Sweep')
        U_sweep{i_k} = U_fluid_o(1:varCol.noDofs,:);
        if ~strcmp(BC,'SHBC')
            error('not implemented due to noDofs')
        end
    end
    if ~useROM
        calculateErrors
    end
    if strcmp(scatteringCase,'Sweep') && size(k,2) > 1
        fprintf('\nTotal time spent on frequency: %12f', toc(t_freq))  
    end
end
if useROM
    varCol.U_sweep = U_sweep;
    varCol.k = k;
end

%% Compute scattered pressure    
tic
if calculateFarFieldPattern && ~useROM
    if ~runTasksInParallel
        fprintf(['\n%-' num2str(stringShift) 's'], 'Computing far-field pattern ... ')
    end
    getFarFielsPoints

    switch method
        case {'IE','ABC'}
            p = calculateScatteredPressure(varCol, U_fluid_o, v, plotFarField);
        case 'IENSG'
            p = calculateScatteredPressureNonSepGeom(varCol, U_fluid_o, v, plotFarField);
        case 'BA'
            p = calculateScatteredPressureBA(varCol, U_fluid_o, v, 0, plotFarField);
        case 'BEM'
            p = calculateScatteredPressureBEM(varCol, U_fluid_o, v, 0, plotFarField);
        case 'KDT'
            switch coreMethod
                case 'linear_FEM'
                    element = varCol.element(:,[1,2,4,3]);

                    indices = or(or(or(or(or(element(:,1) == element(:,2), element(:,1) == element(:,3)), ...
                                                element(:,1) == element(:,4)), ...
                                                element(:,2) == element(:,3)), ...
                                                element(:,2) == element(:,4)), ...
                                                element(:,3) == element(:,4));
                    indices2 = find(indices);
                    for i = 1:numel(indices2)
                        idx = indices2(i);
                        element(idx,:) = [unique(element(idx,:)), NaN];
                    end

                    tri = [element(~indices,1:3);
                           element(~indices,[1,3,4])
                           element(indices,1:3)];

                    P = varCol.controlPts;
                    p = kirchApprTri(tri,P,v,varCol);
%                     trisurf(tri,P(:,1),P(:,2),P(:,3), 'FaceColor', 1.5*[44 77 32]/255)
%                     axis off
%                     axis equal
%                     camlight
%                     
                    %% Find h_max and store results
                    varCol.h_max = max([norm2(P(tri(:,1),:)-P(tri(:,2),:)); norm2(P(tri(:,1),:)-P(tri(:,3),:))]);
                    varCol.dofs = size(unique(tri));
                    varCol.nepw = lambda(1)./varCol.h_max;
                    varCol.noElems = size(tri,1);
                case 'IGA'
                    varCol.h_max = findMaxElementDiameter(varCol.patches);
                    varCol.nepw = lambda(1)./varCol.h_max;
                    p = calculateScatteredPressureKDT(varCol, v, plotFarField);
                    varCol.dofs = varCol.noDofs;
            end
            varCol.tot_time = NaN;
        case 'MFS'
            p = calculateScatteredPressureMFS(varCol, U_fluid_o, v, plotFarField);
        case 'RT'
            switch scatteringCase
                case 'MS'
                    d_vec = varCol.d_vec;
                    p = zeros(size(d_vec,2),1);
                    for i = 1:size(d_vec,2) %874%
%                     parfor i = 1:size(d_vec,2)
                        varColTemp2 = varCol;
                        varColTemp2.d_vec = d_vec(:,i);
                        tic
                        varColTemp2 = createRays(varColTemp2);
                        fprintf('\nCreating rays in %12f seconds.', toc)
                        tic
                        varColTemp2 = traceRays(varColTemp2);    
                        fprintf('\nTracing rays in %12f seconds.', toc)
                        tic        
                        p(i) = calculateScatteredPressureRT(varColTemp2, v(i,:), plotFarField);
                        fprintf('\nFar field in %12f seconds.', toc)
                    end
                otherwise
                    varCol = createRays(varCol);
                    varCol = traceRays(varCol);            
                    p = calculateScatteredPressureRT(varCol, v, plotFarField);
            end
    end
    if ~runTasksInParallel
        fprintf('using %12f seconds.', toc)
    end
    task.results.p = p;
    task.results.abs_p = abs(p);
    task.results.TS = 20*log10(abs(p/P_inc));
    if analyticSolutionExist
        if plotFarField
            p_ref = varCol.farField(v);
%             p_ref = exactKDT(varCol.k,varCol.P_inc,parms.R_o);
        else
            p_ref = varCol.analytic(v);
        end
        task.results.p_ref = p_ref;
        task.results.abs_p_ref = abs(p_ref);
        task.results.TS_ref = 20*log10(abs(p_ref/P_inc));

        task.results.error_pAbs = 100*abs(abs(p_ref)-abs(p))./abs(p_ref);
        task.results.error_p = 100*abs(p_ref-p)./abs(p_ref);
    end
end

if strcmp(method,'KDT') || strcmp(method,'RT')
    return
end
% task.results.U_fluid_o = U_fluid_o;
% if useSolidDomain
%     task.results.varCol_solid = varCol_solid;
% %     task.results.U_solid = U_solid;
% end
% if useInnerFluidDomain
%     task.results.varCol_fluid_i = varCol_fluid_i;
% %     task.results.U_fluid_i = U_fluid_i;
% end
if strcmp(scatteringCase,'Ray')
    plotSolutionAlongRay
end

%% Calculate errors (if analyticSolutionExist) and plot result in Paraview
if ~useROM
    switch scatteringCase
        case {'BI', 'Sweep','Ray'}    
    %         keyboard
            plotError = analyticSolutionExist && ~plotTimeOscillation; 

            tic
            plotResultsInMatlab = 0;
            if plotResultsInMatlab && analyticSolutionExist && strcmp(formulation,'GCBIE')
                plotGalerkinResidual(varCol,U);
                axis equal
                axis off
    %             set(gca, 'Color', 'none');
                title('Error in Galerkin solution')
                colorbar
    %             n_xi = fluid.number(1);
    %             n_eta = fluid.number(2);
    %             [cg_xi, grev_xi] = CauchyGalerkin(fluid.degree(1), n_xi, fluid.knots{1});
    %             [cg_eta, grev_eta] = CauchyGalerkin(fluid.degree(2), n_eta, fluid.knots{2});
    %             % keyboard
    %             n_cp = varCol.noDofs - length(dofsToRemove);
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
            if plotResultsInParaview
                tic
                if ~runTasksInParallel
                    fprintf(['\n%-' num2str(stringShift) 's'], 'Post-processing ... ')
                end
                resultsFolderNameParaview = [subFolderName '/paraviewResults'];
                if ~exist(resultsFolderNameParaview, 'dir')
                    mkdir(resultsFolderNameParaview);
                end          
                vtfFileName = [resultsFolderNameParaview '/' saveName];

    %             noUniqueXiKnots = length(unique(varCol.nurbs.knots{1}));
    %             noUniqueEtaKnots = length(unique(varCol.nurbs.knots{2}));

                fact = 40;

    %             extraXiPts = floor(fact/(varCol.L_gamma/h_max)/2^(M-1)); % .. per element
    %             extraEtaPts = extraXiPts; % .. per element
    %             extraZetaPts = extraXiPts; % .. per element

                extraXiPts = 10; % .. per element
                extraEtaPts = 10; % .. per element
                extraZetaPts = 10; % .. per element
                testFun = @(v) -analytic(v) - P_inc(v);
                varCol.testFun = testFun;
                if boundaryMethod
                    if strcmp(model,'SS_P')
                        options = struct('name',vtfFileName, 'celltype', 'VTK_QUAD', ...
                                    'plotTimeOscillation', plotTimeOscillation, 'plotError', plotError, 'plotTotField', 1, ...
                                    'plotAnalytic', analyticSolutionExist);
                    else
                        options = struct('name',vtfFileName, 'celltype', 'VTK_QUAD', 'plotTimeOscillation', plotTimeOscillation, ...
                                        'plotScalarField',1, 'plotError', plotError, 'plotP_inc', 1, 'plotTotField', 1, ...
                                        'plotTotFieldAbs', 1, ...
                                        'plotErrorFunc', 0, 'plotTestFun', 0, 'plotTestField', 1, 'plotAnalytic', analyticSolutionExist);
                        if strcmp(varCol.applyLoad,'radialPulsation')
                            options.plotP_inc = false;
                            options.plotTotField = false;
                            options.plotTotFieldAbs = false;
                        end
                    end
                elseif strcmp(method, 'IE')
                    options = struct('name',vtfFileName, 'celltype', 'VTK_HEXAHEDRON', 'plotTimeOscillation', plotTimeOscillation, ...
                        'plotScalarField',1, 'plotDisplacementVectors', 1, 'plotError', plotError, 'plotErrorEnergy', 1,...
                        'plotErrorGrad', plotError, 'plotTestFun', 0', 'plotP_inc', 1, 'plotTotField', 1, 'plotTotFieldAbs', 1, 'plotAnalytic', analyticSolutionExist);

                    % plot artificial boundary
                    plotModelInParaview(extractOuterSurface(fluid), extraXiPts, extraEtaPts, extraZetaPts, options, 0, 'artificialBoundary')
                    if strcmp(BC,'SHBC')
                        plotModelInParaview(extractInnerSurface(fluid), extraXiPts, extraEtaPts, extraZetaPts, options, 1, 'scatterer')
                    end
                end
                varCol.isOuterDomain = true;

                if ~boundaryMethod
                    createParaviewFiles(varCol, U_fluid_o, extraXiPts, extraEtaPts, extraZetaPts, options, plotMesh, e3Dss_options);
                end

                if plotMesh
                    createVTKmeshFiles(varCol, U_fluid_o, extraXiPts, extraEtaPts, extraZetaPts, options)
                end

                delta = 0.5;
                xb = [-65,20]+pi*1e-6;
                yb = [-15,15]+pi*1e-6;
                zb = [-10,10]+pi*1e-6;
                options = struct('name',vtfFileName, 'celltype', 'VTK_QUAD', 'plotTimeOscillation', plotTimeOscillation, ...
                                'plotScalarField',1, 'plotError', plotError, 'plotP_inc', 1, 'plotTotField', 1, ...
                                'plotTotFieldAbs', 1, ...
                                'plotErrorFunc', 0, 'plotTestFun', 0, 'plotTestField', 0, 'plotAnalytic', analyticSolutionExist);
                plotTriangulation(varCol,U_fluid_o,delta,xb,yb,zb, options, '_exterior')


                if boundaryMethod
                    fluid2 = fluid;
                    varCol2 = varCol;
                    U_fluid_o2 = U_fluid_o;
                else
                    [fluid2, varCol2] = extractOuterSurface(fluid, varCol);
                    U_fluid_o2 = U_fluid_o(end-task.N*varCol2.noCtrlPts+1:end);
                    varCol2.varColFull = varCol;
                    varCol2.U_full = U_fluid_o;
                end
                if 0
                    R = 50;
                    stringShift = 40;
                    fprintf(['\n%-' num2str(stringShift) 's'], '    Compute solution on sphere ... ')
                    plotOnSphere(varCol2, U_fluid_o2, options, delta, R, zb, '_spherePlot')
                    fprintf('using %12f seconds.', toc)
                end

    %             theta = linspace(0,pi/3,2);

    %             setBCParameters
    %             offset = pi;
    %             Lx = a+offset - (-L-g2-g3-offset);
    %             x = linspace(-L-g2-g3-offset,a+offset, ceil(Lx/delta)+1);
    %             Ly = b+offset - (-b-offset);
    %             y = linspace(-b-offset,b+offset, ceil(Ly/delta)+1);
    %             Lz = c+ht+offset - (-b-offset);
    %             z = linspace(-b-offset,c+ht+offset, ceil(Lz/delta)+1);





    %             if strcmp(model,'BC') || strcmp(model,'BC_P') || strcmp(model,'MS') || strcmp(model,'MS_P') || strcmp(model,'S1') || strcmp(model,'S1_P')
    %                 if boundaryMethod
    %                     fluid2 = fluid;
    %                     varCol2 = varCol;
    %                     U_fluid_o2 = U_fluid_o;
    %                 else
    %                     [fluid2, varCol2] = extractOuterSurface(fluid, varCol);
    %                     U_fluid_o2 = U_fluid_o(end-task.N*varCol2.noCtrlPts+1:end);
    %                     varCol2.varColFull = varCol;
    %                     varCol2.U_full = U_fluid_o;
    %                 end
    %                 if strcmp(model,'BC') || strcmp(model,'BC_P')
    %                     fun = @(xi) objFun2(fluid2,xi,0.8);
    %                     options2 = optimset('TolX', 1e-13);
    %                     xi = fminsearchbnd(fun,0.07,0,0.5,options2);
    %                     xb = [-65,20];
    %                     yb = [-15,15];
    %                     zb = [-10,10];
    %                     xi_arr = [xi, 1-xi];
    %                     delta = 0.11;
    % %                     delta = 1;
    %                 elseif strcmp(model,'MS') || strcmp(model,'MS_P')
    %                     xi_arr = [0, 0.5];
    %                     xb = [-4*R_o-parms.L,4*R_o];
    %                     yb = 3*[-R_o,R_o];
    %                     zb = 3*[-R_o,R_o];
    %                     delta = 0.03;
    %                 elseif strcmp(model,'S1') || strcmp(model,'S1_P')
    %                     xi_arr = [0.25, 0.75];
    %                     xb = [-2*R_o,2*R_o];
    %                     yb = 2*[-R_o,R_o];
    %                     zb = 2*[-R_o,R_o];
    %                     delta = 0.3;
    %                 end
    %                 tic
    %                 plotCrossSection(varCol2, U_fluid_o2, options, [0, 0.8], 1, delta, 1, yb, zb, boundaryMethod, '_crossection1')
    %                 toc
    %                 tic
    %                 plotCrossSection(varCol2, U_fluid_o2, options, [0, 0.7], 1, delta, 1, yb, zb, boundaryMethod, '_crossection2')
    %                 toc
    %                 tic
    %                 options.plotErrorGradient = false;
    %                 options.plotErrorGrad = false;
    %                 options.plotErrorEnergy = false;
    %                 plotCrossSection(varCol2, U_fluid_o2, options, xi_arr, 2, delta, xb, yb, 1, boundaryMethod, '_crossection3')
    %                 toc
    %                 tic
    %                 plotCrossSection(varCol2, U_fluid_o2, options, [0, 0.5], 2, delta, xb, 1, zb, boundaryMethod,'_crossection4')
    %                 toc
    %                 tic
    % %                 
    % %                 
    %             end

    %             if useSolidDomain 
    %                 vtfFileName = [resultsFolderNameParaview '/' saveName '_solid'];
    % 
    %                 noUniqueXiKnots = length(unique(varCol_solid.nurbs.knots{1}));
    %                 noUniqueEtaKnots = length(unique(varCol_solid.nurbs.knots{2}));
    %                 noUniqueZetaKnots = length(unique(varCol_solid.nurbs.knots{3}));
    %                 
    %                 extraZetaPts = 1; % .. per element
    % %                 extraXiPts = 0; % .. per element
    % %                 extraEtaPts = 0; % .. per element
    % %                 extraZetaPts = 0; % .. per element
    %                 options = struct('name',vtfFileName, 'celltype', 'VTK_HEXAHEDRON', 'plotTimeOscillation', plotTimeOscillation, ...
    %                                 'plotErrorGrad', 1, 'plotDisplacementVectors',1, 'plotErrorEnergy', 1, ...                   
    %                 'plotSphericalStress_rr',1, 'plotError', 1,'plotVonMisesStress',1, ...
    %                  'plotStressXX',1,...
    %                  'plotStressYY',1,...
    %                  'plotStressZZ',1,...
    %                  'plotStressYZ',1,...
    %                  'plotStressXZ',1,...
    %                  'plotStressXY',1); 
    % 
    %                 varCol_solid.isOuterDomain = false;
    %                 createParaviewFiles(varCol_solid, U_solid, extraXiPts, extraEtaPts, extraZetaPts, options, plotMesh, e3Dss_options);
    %             end
    % 
    %             if useInnerFluidDomain 
    %                 vtfFileName = [resultsFolderNameParaview '/' saveName '_inner'];
    % 
    %                 noUniqueXiKnots = length(unique(varCol_fluid_i.nurbs.knots{1}));
    %                 noUniqueEtaKnots = length(unique(varCol_fluid_i.nurbs.knots{2}));
    %                 noUniqueZetaKnots = length(unique(varCol_fluid_i.nurbs.knots{3}));
    %                 
    %                 extraZetaPts = floor(fact/(R_i(1)*2/h_max)/2^(M-1)); % .. per element
    % %                 extraZetaPts = floor(fact/(2*pi*b)/(noUniqueZetaKnots-1)); % .. per element
    % %                 extraXiPts = 0; % .. per element
    % %                 extraEtaPts = 0; % .. per element
    % %                 extraZetaPts = 0; % .. per element
    %                 options = struct('name',vtfFileName, 'celltype', 'VTK_HEXAHEDRON', 'plotTimeOscillation', plotTimeOscillation, ...
    %                         'plotTotField', 1, 'plotErrorGrad', 1, 'plotDisplacementVectors', 1, 'plotErrorEnergy', 1, ...
    %                         'plotGradientVectors', 0, 'plotError', 1, 'plotTotFieldAbs', ~plotTimeOscillation);
    % 
    %                 gP_inc = @(v) zeros(length(v),3);
    % 
    %                 varCol_fluid_i.isOuterDomain = false;
    %                 createParaviewFiles(varCol_fluid_i, U_fluid_i, extraXiPts, extraEtaPts, extraZetaPts, options, plotMesh, e3Dss_options);
    %             end
            end
    end
end

if storeFullVarCol
    varColTemp = varCol;
else
    varColTemp.alpha = varCol.alpha;
    varColTemp.beta = varCol.beta;
    varColTemp.k = k;
    varColTemp.f = k*varCol.c_f(1)/(2*pi);
    if ~strcmp(method,'RT')
        varColTemp.dofs = varCol.dofs;
        varColTemp.dofsAlg = (varCol.dofs)^(1/3);
        varColTemp.nepw = varCol.nepw;
        varColTemp.tot_time = varCol.tot_time;
        varColTemp.noElems = varCol.noElems;
        if storeSolution
            varColTemp.U = varCol.U;
        end
        if useROM
            varColTemp.U_sweep = U_sweep;
        end
        if isfield(varCol,'N')
            varColTemp.N = varCol.N;
        end
        varColTemp.h_max = varCol.h_max;
        if isfield(varCol,'timeBuildSystem')
            varColTemp.timeBuildSystem = varCol.timeBuildSystem;
        end
        if isfield(varCol,'timeSolveSystem')
            varColTemp.timeSolveSystem = varCol.timeSolveSystem;
        end
    end
    varColTemp.analyticSolutionExist = varCol.analyticSolutionExist;
end
task.varCol = varColTemp;

fprintf('\nTotal time spent on task: %12f', toc(t_start))  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  
% cond_number = task.results.cond_number;
% p = task.results.p;
% abs_p = task.results.abs_p;
% TS = task.results.TS;
% 
% p_ref = task.results.p_ref;
% abs_p_ref = task.results.abs_p_ref;
% TS_ref = task.results.TS_ref;
% 
% error_pAbs = task.results.error_pAbs;
% error_p = task.results.error_p;
% energyError = task.results.energyError;
% L2Error = task.results.L2Error;
% H1Error = task.results.H1Error;
% % surfaceError = task.results.surfaceError;
% k = task.varCol.k;
% f = task.varCol.f;
% dofs = task.varCol.dofs;
% nepw = task.varCol.nepw;
% tot_time = task.varCol.tot_time;
% N = task.varCol.N;
% h_max = task.varCol.h_max;
% timeBuildSystem = task.varCol.timeBuildSystem;
% analyticSolutionExist = task.varCol.analyticSolutionExist;
