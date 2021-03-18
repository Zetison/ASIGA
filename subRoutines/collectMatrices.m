function [task,FF,A0,A1,A2,A4] = collectMatrices(task)
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
A0 = sparse(noRows_tot,noCols_tot);
A1 = sparse(noRows_tot,noCols_tot);
A2 = sparse(noRows_tot,noCols_tot);
A4 = sparse(noRows_tot,noCols_tot);
FF = zeros(noRows_tot,size(task.varCol{1}.FF,2));
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
            if strcmp(task.misc.method,'BA')
                massMatrixScale = eqScale;
            else
                massMatrixScale = -eqScale/c_f^2;
            end
        case 'solid'
            eqScale = 1;
            if strcmp(task.misc.method,'BA')
                massMatrixScale = eqScale;
            else
                massMatrixScale = -eqScale*rho;
            end
    end
    if isfield(task.varCol{i},'A_K')
        A0(Aindices{i,1},Aindices{i,2}) = eqScale*task.varCol{i}.A_K; 
    end
    
    if isfield(task.varCol{i},'A_M')
        A2(Aindices{i,1},Aindices{i,2}) = A2(Aindices{i,1},Aindices{i,2}) + massMatrixScale*task.varCol{i}.A_M; 
    end
    if isfield(task.varCol{i},'A_C') && i > 1
        Cindices1 = Aindices{i,1}(end)+(1:size(task.varCol{i}.A_C,1));
        Cindices2 = Aindices{i,2}(1)-1+(1:size(task.varCol{i}.A_C,2));
        switch task.varCol{i-1}.media
            case 'fluid'
                A2(Cindices1,Cindices2) = A2(Cindices1,Cindices2) + task.varCol{i}.A_C;
            case 'solid'
                A0(Cindices1,Cindices2) = A0(Cindices1,Cindices2) + task.varCol{i}.A_C;
        end
        switch task.varCol{i}.media
            case 'fluid'
                A2(Cindices2,Cindices1) = A2(Cindices2,Cindices1) + task.varCol{i}.A_C.'; 
            case 'solid'
                A0(Cindices2,Cindices1) = A0(Cindices2,Cindices1) + task.varCol{i}.A_C.'; 
        end
    end
    if isfield(task.varCol{i},'FF')
        FF(Aindices{i,1},:) = eqScale*task.varCol{i}.FF;
    end
end
if strcmp(task.misc.method,'IE') || strcmp(task.misc.method,'IENSG')
    if isfield(task.varCol{1},'Ainf')
    	A0(AindicesInf,AindicesInf) = A0(AindicesInf,AindicesInf) + task.varCol{1}.Ainf/task.varCol{1}.rho; 
        if task.rom.useROM
            A1(AindicesInf,AindicesInf) = A1(AindicesInf,AindicesInf) + task.varCol{1}.Ainf1/task.varCol{1}.rho/task.varCol{1}.c_f; 
            A2(AindicesInf,AindicesInf) = A2(AindicesInf,AindicesInf) + task.varCol{1}.Ainf2/task.varCol{1}.rho/task.varCol{1}.c_f^2; 
        end  
    end
end

if ~(strcmp(task.misc.method,'BEM') && strcmp(task.misc.formulation(1),'C'))
    A0(allDofsToRemove,:) = [];
    A1(allDofsToRemove,:) = [];
    A2(allDofsToRemove,:) = [];
    A4(allDofsToRemove,:) = [];
    FF(allDofsToRemove,:) = [];
end
A0(:,allDofsToRemove) = [];
A1(:,allDofsToRemove) = [];
A2(:,allDofsToRemove) = [];
A4(:,allDofsToRemove) = [];
