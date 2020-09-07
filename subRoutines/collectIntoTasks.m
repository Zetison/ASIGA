variables = whos;
protectedVariables = {'variables','i','j','counter','studyName','studies','studiesCol','loopParameters'};
for i = 1:length(variables)
    if ~any(ismember(variables(i).name,protectedVariables))
        task.(variables(i).name) = eval(variables(i).name);
    end
end
if task.useROM
    task.storeSolution = true;
end
if (~isnan(alpha_s(1)) || ~isnan(beta_s(1))) && (strcmp(scatteringCase,'MS') || strcmp(scatteringCase,'Sweep'))
    error(['For monostatic scattering alpha_s and beta_s should not be given (they should be defined through alpha and beta). ' ...
           'Note that alpha_s = alpha and beta_s = beta.'])
end
if (isnan(alpha_s(1)) || isnan(beta_s(1))) && strcmp(scatteringCase,'BI') && strcmp(applyLoad,'planeWave')
    error('Incident direction is not set: alpha_s = NaN and/or beta_s = NaN')
end
if solveForPtot && ~(strcmp(method,'BEM') || strcmp(method,'BA'))
    error('solveForPtot should can only be used with method = BEM or method = BA')
end
loopParametersArr = cell(length(loopParameters),1);

taskNames = fieldnames(task);
for j = 1:numel(taskNames)
    idx = find(strcmp(taskNames{j},loopParameters));
    if ~isempty(idx)
        loopParametersArr{idx} = task.(taskNames{j});
        if ischar(loopParametersArr{idx}) % ensure cell type
            loopParametersArr{idx} = loopParametersArr(idx);
        end
    elseif iscell(task.(taskNames{j})) && ~strcmp(taskNames{j},'varCol') % remove redundant cell type
        temp = task.(taskNames{j});
        task.(taskNames{j}) = temp{1};
    end
end

studies(counter).loopParameters = loopParameters;
studies(counter).loopParametersArr = loopParametersArr;
studies(counter).runTasksInParallel = runTasksInParallel;
if exist('basisROMcell','var')
    studies(counter).basisROMcell = basisROMcell;
    studies(counter).k_ROM = k_ROM;
    studies(counter).noVecsArr = noVecsArr;
end

studies(counter).tasks = createTasks([], 1, task, 1, loopParameters, loopParametersArr);
studies(counter).postPlot = postPlot;
if isempty(subFolderName)
    subFolderName = model;
end
studies(counter).subFolderName = subFolderName;

counter = counter + 1;

