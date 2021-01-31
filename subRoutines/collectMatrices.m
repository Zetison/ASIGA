function [varCol,FF,A0,A1,A2,A4] = collectMatrices(varCol,task)

noDomains = numel(varCol);
Aindices = cell(noDomains,2);
noCols_tot = 0;
noRows_tot = 0;
if isfield(varCol{1},'A_K')
    allDofsToRemove = [];
    for i = noDomains:-1:1
        Aindices{i,1} = noRows_tot+(1:size(varCol{i}.A_K,1));
        Aindices{i,2} = noCols_tot+(1:size(varCol{i}.A_K,2));
        allDofsToRemove = [allDofsToRemove (varCol{i}.dofsToRemove+noCols_tot)];
        noRows_tot = noRows_tot + size(varCol{i}.A_K,1);
        noCols_tot = noCols_tot + size(varCol{i}.A_K,2);
    end
    if strcmp(task.method,'IE') || strcmp(task.method,'IENSG')
        AindicesInf = noCols_tot - varCol{1}.noDofs+(1:varCol{1}.noDofs_new);
        noCols_tot = noCols_tot - varCol{1}.noDofs + varCol{1}.noDofs_new;
        noRows_tot = noRows_tot - varCol{1}.noDofs + varCol{1}.noDofs_new;
    end
    varCol{1}.noCols_tot = noCols_tot;
    varCol{1}.noRows_tot = noRows_tot;
    varCol{1}.Aindices = Aindices;
    varCol{1}.allDofsToRemove = allDofsToRemove;
else
    noCols_tot = varCol{1}.noCols_tot;
    noRows_tot = varCol{1}.noRows_tot;
    Aindices = varCol{1}.Aindices;
    allDofsToRemove = varCol{1}.allDofsToRemove;
end
% Collect all matrices
A0 = sparse(noRows_tot,noCols_tot);
A1 = sparse(noRows_tot,noCols_tot);
A2 = sparse(noRows_tot,noCols_tot);
A4 = sparse(noRows_tot,noCols_tot);
FF = zeros(noRows_tot,varCol{1}.noRHSs);
for i = 1:noDomains
    switch varCol{i}.media
        case 'fluid'
            eqScale = 1/varCol{i}.rho;
            massMatrixScale = -eqScale/varCol{i}.c_f^2;
        case 'solid'
            eqScale = 1;
            massMatrixScale = -eqScale*varCol{i}.rho;
    end
    if isfield(varCol{i},'A_K')
        A0(Aindices{i,1},Aindices{i,2}) = eqScale*varCol{i}.A_K; 
    end
    
    if isfield(varCol{i},'A_M')
        A2(Aindices{i,1},Aindices{i,2}) = A2(Aindices{i,1},Aindices{i,2}) + massMatrixScale*varCol{i}.A_M; 
    end
    if isfield(varCol{i},'A_C') && i > 1
        Cindices1 = Aindices{i,1}(end)+(1:size(varCol{i}.A_C,1));
        Cindices2 = Aindices{i,2}(1)-1+(1:size(varCol{i}.A_C,2));
        switch varCol{i-1}.media
            case 'fluid'
                A2(Cindices1,Cindices2) = A2(Cindices1,Cindices2) + varCol{i}.A_C;
            case 'solid'
                A0(Cindices1,Cindices2) = A0(Cindices1,Cindices2) + varCol{i}.A_C;
        end
        switch varCol{i}.media
            case 'fluid'
                A2(Cindices2,Cindices1) = A2(Cindices2,Cindices1) + varCol{i}.A_C.'; 
            case 'solid'
                A0(Cindices2,Cindices1) = A0(Cindices2,Cindices1) + varCol{i}.A_C.'; 
        end
    end
    if isfield(varCol{i},'FF')
        FF(Aindices{i,1},:) = eqScale*varCol{i}.FF;
    end
end
if strcmp(task.method,'IE') || strcmp(task.method,'IENSG')
    if isfield(varCol{i},'Ainf')
    	A0(AindicesInf,AindicesInf) = A0(AindicesInf,AindicesInf) + varCol{1}.Ainf/varCol{1}.rho; 
        if task.useROM
            A1(AindicesInf,AindicesInf) = A1(AindicesInf,AindicesInf) + varCol{1}.Ainf1/varCol{1}.rho/varCol{i}.c_f; 
            A2(AindicesInf,AindicesInf) = A2(AindicesInf,AindicesInf) + varCol{1}.Ainf2/varCol{1}.rho/varCol{i}.c_f^2; 
        end  
    end
end

if ~(strcmp(task.method,'BEM') && strcmp(task.formulation(1),'C'))
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
