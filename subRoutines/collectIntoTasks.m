
structs = {'varCol','misc','msh','prePlot','postPlot','sol','err','ffp','para','iem','pml','mfs','rom'};
for i = 1:numel(structs)
    task.(structs{i}) = eval(structs{i});
end
if task.rom.useROM
    task.misc.storeSolution = true;
end
if (~isnan(ffp.alpha_s(1)) || ~isnan(ffp.beta_s(1))) && (strcmp(misc.scatteringCase,'MS'))
    error(['For monostatic scattering alpha_s and beta_s should not be given (they should be defined through alpha and beta). ' ...
           'Note that alpha_s = alpha and beta_s = beta.'])
end
if (isnan(ffp.alpha_s(1)) || isnan(ffp.beta_s(1))) && strcmp(misc.scatteringCase,'BI') && strcmp(misc.applyLoad,'planeWave')
    error('Incident direction is not set: alpha_s = NaN and/or beta_s = NaN')
end
if misc.solveForPtot && ~(strcmp(misc.method,'BEM') || strcmp(misc.method,'BA'))
    error('solveForPtot should can only be used with method = BEM or method = BA')
end
loopParametersArr = cell(length(loopParameters),1);

taskNames = fieldnames(task);
allTaskNames = cell(0,0);
for j = 1:numel(taskNames)
    taskName_j = task.(taskNames{j});
    if isstruct(taskName_j)
        subTaskNames = fieldnames(taskName_j);
        for i = 1:numel(subTaskNames)
            allTaskNames{end+1,1} = [taskNames{j} '.' subTaskNames{i}];
        end
    else
        allTaskNames(end+1,1) = taskNames(j);
    end
end
for j = 1:numel(allTaskNames)
    idx = find(strcmp(allTaskNames{j},loopParameters));
    if ~isempty(idx)
        loopParametersArr{idx} = eval(['task.' allTaskNames{j}]);
        if ischar(loopParametersArr{idx}) % ensure cell type
            loopParametersArr{idx} = loopParametersArr(idx);
        end
    elseif iscell(eval(['task.' allTaskNames{j}])) && isempty(strfind(allTaskNames{j},'postPlot')) && isempty(strfind(allTaskNames{j},'prePlot')) ...
            && ~strcmp(allTaskNames{j},'varCol') % remove redundant cell type
        temp = eval(['task.' allTaskNames{j}]);
        eval(['task.' allTaskNames{j} ' = temp{1};'])
    end
end

studies(counter).loopParameters = loopParameters;
studies(counter).loopParametersArr = loopParametersArr;
studies(counter).runTasksInParallel = runTasksInParallel;
if exist('basisROMcell','var')
    studies(counter).basisROMcell = basisROMcell;
    studies(counter).omega_ROM = omega_ROM;
    studies(counter).noVecsArr = noVecsArr;
end

studies(counter).tasks = createTasks([], 1, task, 1, loopParameters, loopParametersArr);
studies(counter).postPlot = postPlot;
if isempty(subFolderName)
    subFolderName = misc.model;
end
studies(counter).subFolderName = subFolderName;

counter = counter + 1;

