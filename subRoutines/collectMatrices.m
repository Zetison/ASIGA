function task = collectMatrices(task,collectLHS,collectRHS,buildA)
if nargin < 2
    collectLHS = true;
end
if nargin < 3
    collectRHS = true;
end
if nargin < 4
    buildA = true;
end
if isfield(task,'sol')
    useCSLP = strcmp(task.sol.preconditioner,'CSLP');
else
    useCSLP = false;
end
switch task.misc.method
    case {'IE','IENSG','BEM','BA','ABC','MFS','PML'}
        tic
        if task.misc.printLog
            fprintf(['\n%-' num2str(task.misc.stringShift) 's'], 'Collecting matrices ... ')
        end
        useA = ~(task.noDomains == 1 && strcmp(task.misc.method,'BEM'));
        if useA
            noDomains = numel(task.varCol);
            Aindices = cell(noDomains,2);
            noCols_tot = 0;
            noRows_tot = 0;
            if isfield(task.varCol{1},'A_K') || isfield(task.varCol{1},'A_M') || strcmp(task.misc.method,'IENSG')
                allDofsToRemove = [];
                for i = noDomains:-1:1
                    if isfield(task.varCol{i},'A_K')
                        noRows = size(task.varCol{i}.A_K,1);
                        noCols = size(task.varCol{i}.A_K,2);
                    elseif isfield(task.varCol{i},'A_M')
                        noRows = size(task.varCol{i}.A_M,1);
                        noCols = size(task.varCol{i}.A_M,2);
                    elseif strcmp(task.misc.method,'IENSG')
                        noRows = task.varCol{i}.noDofs;
                        noCols = task.varCol{i}.noDofs;
                    else
                        error('No matrix found for this domain')
                    end
                    Aindices{i,1} = noRows_tot+(1:noRows);
                    Aindices{i,2} = noCols_tot+(1:noCols);
                    allDofsToRemove = [allDofsToRemove (task.varCol{i}.dofsToRemove+noCols_tot)];
                    noRows_tot = noRows_tot + noRows;
                    noCols_tot = noCols_tot + noCols;
                end
                if strcmp(task.misc.method,'IE') || strcmp(task.misc.method,'IENSG')
                    AindicesInf = noCols_tot - task.varCol{1}.noDofs+(1:task.varCol{1}.noDofs_new);
                    noCols_tot = noCols_tot - task.varCol{1}.noDofs + task.varCol{1}.noDofs_new;
                    noRows_tot = noRows_tot - task.varCol{1}.noDofs + task.varCol{1}.noDofs_new;
                elseif strcmp(task.misc.method,'ABC') 
                    AindicesInf = noCols_tot - task.varCol{1}.noDofs + (1:task.varCol{1}.noDofs);
                end
                task.varCol{1}.noCols_tot = noCols_tot;
                task.varCol{1}.noRows_tot = noRows_tot;
                task.varCol{1}.Aindices = Aindices;
                task.varCol{1}.allDofsToRemove = allDofsToRemove;
            else
                noCols_tot = task.varCol{1}.noCols_tot;
                noRows_tot = task.varCol{1}.noRows_tot;
                Aindices = task.varCol{1}.Aindices;
                allDofsToRemove = task.varCol{1}.allDofsToRemove;
            end
            % Collect all matrices
            if collectLHS
                task.A0 = sparse(noRows_tot,noCols_tot);
                task.A1 = sparse(noRows_tot,noCols_tot);
                task.A2 = sparse(noRows_tot,noCols_tot);
                task.A4 = sparse(noRows_tot,noCols_tot);
                task.P_right = speye(noRows_tot,noCols_tot);
                task.P_rightinv = speye(noRows_tot,noCols_tot);
            end
            if collectRHS
                task.FF = complex(zeros(noRows_tot,size(task.varCol{1}.FF,2)));
            end
            for i = 1:noDomains
                if isfield(task.varCol{i},'rho')
                    rho = task.varCol{i}.rho;
                else
                    rho = 1;
                end
                switch task.varCol{i}.media
                    case 'fluid'
                        if isfield(task.varCol{i},'c_f')
                            c_f = task.varCol{i}.c_f;
                        else
                            c_f = 1;
                        end
                        eqScale = 1/rho;
                        if strcmp(task.misc.method,'BA') || strcmp(task.misc.method,'interp')
                            massMatrixScale = eqScale;
                        else
                            massMatrixScale = -eqScale/c_f^2;
                        end
                    case 'solid'
                        eqScale = 1;
                        if strcmp(task.misc.method,'BA') || strcmp(task.misc.method,'interp')
                            massMatrixScale = eqScale;
                        else
                            massMatrixScale = -eqScale*rho;
                        end
                end
                if collectLHS
                    if isfield(task.varCol{i},'A_K')
                        if strcmp(task.varCol{i}.media,'solid') && task.misc.symmetric
                            task.A2(Aindices{i,1},Aindices{i,2}) = eqScale*task.varCol{i}.A_K; 
                        else
                            task.A0(Aindices{i,1},Aindices{i,2}) = eqScale*task.varCol{i}.A_K; 
                        end
                        task.varCol{i} = rmfield(task.varCol{i},'A_K');
                    end
                    
                    if isfield(task.varCol{i},'A_M')
                        if task.misc.symmetric
                            switch task.varCol{i}.media
                                case 'fluid'
                                    task.A2(Aindices{i,1},Aindices{i,2}) = task.A2(Aindices{i,1},Aindices{i,2}) + massMatrixScale*task.varCol{i}.A_M; 
                                case 'solid'
                                    task.A4(Aindices{i,1},Aindices{i,2}) = task.A4(Aindices{i,1},Aindices{i,2}) + massMatrixScale*task.varCol{i}.A_M; 
                            end
                        else
                            task.A2(Aindices{i,1},Aindices{i,2}) = task.A2(Aindices{i,1},Aindices{i,2}) + massMatrixScale*task.varCol{i}.A_M; 
                        end
                        if ~useCSLP
                            task.varCol{i} = rmfield(task.varCol{i},'A_M');
                        end
                    end
                    if isfield(task.varCol{i},'A_C') && i > 1
                        Cindices1 = Aindices{i,1}(end)+(1:size(task.varCol{i}.A_C,1));
                        Cindices2 = Aindices{i,2}(1)-1+(1:size(task.varCol{i}.A_C,2));
                        if task.misc.symmetric
                            task.A2(Cindices1,Cindices2) = task.A2(Cindices1,Cindices2) + task.varCol{i}.A_C;
                            task.A2(Cindices2,Cindices1) = task.A2(Cindices2,Cindices1) + task.varCol{i}.A_C.'; 
                        else
                            switch task.varCol{i-1}.media
                                case 'fluid'
                                    task.A2(Cindices1,Cindices2) = task.A2(Cindices1,Cindices2) + task.varCol{i}.A_C;
                                case 'solid'
                                    task.A0(Cindices1,Cindices2) = task.A0(Cindices1,Cindices2) + task.varCol{i}.A_C;
                            end
                            switch task.varCol{i}.media
                                case 'fluid'
                                    task.A2(Cindices2,Cindices1) = task.A2(Cindices2,Cindices1) + task.varCol{i}.A_C.'; 
                                case 'solid'
                                    task.A0(Cindices2,Cindices1) = task.A0(Cindices2,Cindices1) + task.varCol{i}.A_C.'; 
                            end
                        end
                        task.varCol{i} = rmfield(task.varCol{i},'A_C');
                    end
                    if isfield(task.misc,'omega')
                        omega_mean = task.omega_mean;
                        if i > 1 && strcmp(task.varCol{2}.media, 'solid') && strcmp(task.varCol{1}.media, 'fluid')
                            rho = task.varCol{1}.rho;
                            task.P_right(Aindices{i,1},Aindices{i,2}) = omega_mean^2*rho.*speye(task.varCol{i}.noDofs);
                            task.P_rightinv(Aindices{i,1},Aindices{i,2}) = 1/(omega_mean^2*rho).*speye(task.varCol{i}.noDofs);
                        end
                    end
                end
                if collectRHS && isfield(task.varCol{i},'FF')
                    if task.misc.symmetric && strcmp(task.varCol{i}.media,'solid')
                        omega = task.misc.omega;
                        task.FF(Aindices{i,1},:) = omega.^2.*eqScale.*task.varCol{i}.FF;
                    else
                        task.FF(Aindices{i,1},:) = eqScale*task.varCol{i}.FF;
                    end
                    task.varCol{i} = rmfield(task.varCol{i},'FF');
                end
            end
            
            if collectLHS && (strcmp(task.misc.method,'IE') || strcmp(task.misc.method,'IENSG') || strcmp(task.misc.method,'ABC'))
                if isfield(task.varCol{1},'Ainf')
    	            task.A0(AindicesInf,AindicesInf) = task.A0(AindicesInf,AindicesInf) + task.varCol{1}.Ainf/task.varCol{1}.rho; 
                    task.varCol{1} = rmfield(task.varCol{1},'Ainf');
                    if task.rom.useROM
                        task.A1(AindicesInf,AindicesInf) = task.A1(AindicesInf,AindicesInf) + task.varCol{1}.Ainf1/task.varCol{1}.rho/task.varCol{1}.c_f; 
                        task.varCol{1} = rmfield(task.varCol{1},'Ainf1');
                        task.A2(AindicesInf,AindicesInf) = task.A2(AindicesInf,AindicesInf) + task.varCol{1}.Ainf2/task.varCol{1}.rho/task.varCol{1}.c_f^2; 
                        task.varCol{1} = rmfield(task.varCol{1},'Ainf2');
                    end  
                end
            end
            if ~((strcmp(task.misc.method,'BEM') || strcmp(task.misc.method,'interp')) && strcmp(task.misc.formulation(1),'C'))
                if collectLHS
                    task.A0(allDofsToRemove,:) = [];
                    task.A1(allDofsToRemove,:) = [];
                    task.A2(allDofsToRemove,:) = [];
                    task.A4(allDofsToRemove,:) = [];
                    task.P_right(allDofsToRemove,:) = [];
                    task.P_rightinv(allDofsToRemove,:) = [];
                    if useCSLP
                        task.varCol{1}.A_M(allDofsToRemove,:) = [];
                    end
                end
                if collectRHS
                    task.FF(allDofsToRemove,:) = [];
                end
            end
            if collectLHS
                task.A0(:,allDofsToRemove) = [];
                task.A1(:,allDofsToRemove) = [];
                task.A2(:,allDofsToRemove) = [];
                task.A4(:,allDofsToRemove) = [];
                task.P_right(:,allDofsToRemove) = [];
                task.P_rightinv(:,allDofsToRemove) = [];
                if useCSLP
                    task.varCol{1}.A_M(:,allDofsToRemove) = [];
                end
            end
            if buildA
                if strcmp(task.misc.method,'BA')
                    task.A = task.A2;
                else
                    omega = task.misc.omega;
                    task.A = task.A0 + omega*task.A1 + omega^2*task.A2 + omega^4*task.A4;
                end
            end
        else
            if useCSLP
                error('The mass matrix is not available for this method and the CSLP is thus not applicable')
            end
            allDofsToRemove = task.varCol{1}.dofsToRemove;
            if ~(strcmp(task.misc.method,'BEM') && strcmp(task.misc.formulation(1),'C'))
                task.varCol{1}.A_K(allDofsToRemove,:) = [];
                task.varCol{1}.FF(allDofsToRemove,:) = [];
            end
            task.varCol{1}.A_K(:,allDofsToRemove) = [];
        end
        task.timeCollectingMatrices = toc;
        if task.misc.printLog
            fprintf('using %12f seconds.', task.timeCollectingMatrices)
        end
end