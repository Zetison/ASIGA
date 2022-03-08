
structs = {'varCol','misc','msh','prePlot','sol','err','ffp','para','iem','pml','bem','mfs','rt','rom'};
for i = 1:numel(structs)
    task.(structs{i}) = eval(structs{i});
end
if task.rom.useROM
    task.misc.storeSolution = true;
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
childrenParameters = cell(size(connectedParameters));
for i = 1:numel(connectedParameters)
    childrenParameters{i} = cell(size(connectedParameters{i}));
    for j = 1:numel(connectedParameters{i})
        eval(['childrenParameters{i}{j} = task.' connectedParameters{i}{j} ';'])
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
            && isempty(strfind(allTaskNames{j},'para')) && isempty(strfind(allTaskNames{j},'rom')) && ~strcmp(allTaskNames{j},'varCol') % remove redundant cell type
        temp = eval(['task.' allTaskNames{j}]);
        eval(['task.' allTaskNames{j} ' = temp{1};'])
    end
end

if saveStudies && runTasksInParallel
    error('This may generate errors due to large memory consumption of studiesCol')
end
studies(counter).loopParameters = loopParameters;
studies(counter).loopParametersArr = loopParametersArr;
studies(counter).runTasksInParallel = runTasksInParallel;
studies(counter).noCoresToUse = noCoresToUse;
studies(counter).saveStudies = saveStudies;

task.rom.useGDP = false;
for i = 1:numel(rom.basisROMcell)
    if strcmp(rom.basisROMcell{i},'DGP')
        task.rom.useDGP = true;
    end
end
if task.rom.useROM
    studies(counter).basisROMcell = task.rom.basisROMcell;
    studies(counter).omega_ROM = task.rom.omega_ROM;
    studies(counter).noVecsArr = task.rom.noVecsArr;
end

studies(counter).tasks = createTasks([], 1, task, 1, loopParameters, loopParametersArr, connectedParameters, childrenParameters);
if isempty(studies(counter).tasks)
    error('loopParameters contains invalid parameters')
end
studies(counter).postPlot = postPlot;
if isempty(subFolderName)
    subFolderName = misc.model;
end
studies(counter).subFolderName = subFolderName;

counter = counter + 1;

