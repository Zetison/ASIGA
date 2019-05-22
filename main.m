if exist('../availableStudies.m', 'file')
    clear all %#ok
    printAndPlotResults = true;
    availableStudies
else
    printAndPlotResults = false;
end
studies = getTasks(studyName);
startMatlabPool
% Create folders
folderName = '../results';
if ~exist(folderName, 'dir')
    mkdir(folderName);
end

folderName = [folderName '/_studies'];
if ~exist(folderName, 'dir')
    mkdir(folderName);
end
subFolderName = [folderName '/' studyName];
i = 1;
while exist(subFolderName, 'dir')
    i = i + 1;
    subFolderName = [folderName '/' studyName num2str(i)];
end
mkdir(subFolderName);

for study_i = 1:numel(studies)    
    loopParameters = studies(study_i).loopParameters;
    loopParametersArr = studies(study_i).loopParametersArr;
    noTasks = length(studies(study_i).tasks);
    
    tasks = studies(study_i).tasks;
    runTasksInParallel = studies(study_i).runTasksInParallel;
    if runTasksInParallel
%         for i_task = 1:length(tasks)
        parfor i_task = 1:length(tasks)
            tasks(i_task).task = main_sub(tasks(i_task).task,loopParameters,runTasksInParallel,subFolderName);
            fprintf('\nCompleted task %d/%d in study %d/%d\n\n', i_task, noTasks, study_i,length(studies)) 
        end
        studies(study_i).tasks = tasks;
        studies(study_i).subFolderName = subFolderName;
        save([subFolderName '/_studies'], 'studies', '-v7.3')
    else
        for i_task = 1:noTasks
            tasks(i_task).task = main_sub(tasks(i_task).task,loopParameters,runTasksInParallel,subFolderName);
            task = tasks(i_task).task;
            if task.plot3Dgeometry || task.plot2Dgeometry
                return
            end
            if task.useROM
                computeROMsolution
            end
            studies(study_i).tasks(i_task,1).task = tasks(i_task,1).task;
            studies(study_i).subFolderName = subFolderName;
            fprintf('\nCompleted task %d/%d in study %d/%d\n\n', i_task, noTasks, study_i,length(studies)) 
            save([subFolderName '/_studies'], 'studies')
        end
        if task.useROM
            studies(study_i).loopParameters{end+1} = 'noVecs';
            studies(study_i).loopParametersArr{end+1} = noVecsArr;
            studies(study_i).loopParameters{end+1} = 'basisROM';
            studies(study_i).loopParametersArr{end+1} = basisROMcell;
            studies(study_i).tasks = tasks(:);
        end
        save([subFolderName '/_studies'], 'studies')
    end
end

plotFileName = ['plotResults_' studyName];
if printAndPlotResults && exist(['postProcessing/' plotFileName], 'file')
    close all
    eval(plotFileName)
    if isfield(options,'subFolderName') && ~isempty(options.subFolderName)
        save([options.subFolderName, '/_studies'],'studies')
    end
end

