%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Main file of ASIGA
% Author: Jon Vegard Venås
% E-mail: JonVegard.Venas@sintef.no
% Institute: SINTEF Digital
% Release: 2
% Release date: 29/07/2020

%% Extract tasks
if ~exist('../e3Dss', 'dir')
    error('The e3Dss repository must be downloaded (can be obtained from GitHub: https://github.com/Zetison/e3Dss) and placed at the same level as ASIGA (as ../e3Dss)');
end

startup

if exist('../availableStudies.m', 'file')
    clear all %#ok
    availableStudies
end
if ~iscell(studyName)
    studyName = {studyName};
end
studies = getTasks(studyName);
startMatlabPool

%% Create folders
folderName = '../results';
if ~exist(folderName, 'dir')
    mkdir(folderName);
end

folderName = [folderName '/_studies'];
if ~exist(folderName, 'dir')
    mkdir(folderName);
end
subFolderName = [folderName '/' studyName{1}];
i = 1;
while exist(subFolderName, 'dir')
    i = i + 1;
    subFolderName = [folderName '/' studyName{1} '_' num2str(i)];
end
mkdir(subFolderName);

%% Perform studies
t_start_study = tic;
for i_study = 1:numel(studies)    
    loopParameters = studies(i_study).loopParameters;
    loopParametersArr = studies(i_study).loopParametersArr;
    noTasks = length(studies(i_study).tasks);
    
    tasks = studies(i_study).tasks;
    runTasksInParallel = studies(i_study).runTasksInParallel;
    if runTasksInParallel
%         for i_task = 1:length(tasks)
        parfor i_task = 1:length(tasks)
            tasks(i_task).task = main_sub(tasks(i_task).task,loopParameters,runTasksInParallel,subFolderName);
            fprintf('\nStudy %s: Completed task %d/%d in study %d/%d\n\n', studyName{1}, i_task, noTasks, i_study,length(studies)) 
        end
        studies(i_study).tasks = tasks;
        studies(i_study).subFolderName = subFolderName;
        save([subFolderName '/_studies'], 'studies', '-v7.3')
    else
        for i_task = 1:noTasks
            tasks(i_task).task = main_sub(tasks(i_task).task,loopParameters,runTasksInParallel,subFolderName);
            task = tasks(i_task).task;
            if (task.prePlot.plot3Dgeometry || task.prePlot.plot2Dgeometry) && task.prePlot.abortAfterPlotting
                return
            end
            if task.useROM
                computeROMsolution
            end
            studies(i_study).tasks(i_task,1).task = tasks(i_task,1).task;
            studies(i_study).subFolderName = subFolderName;
            fprintf('\nStudy %s: Completed task %d/%d in study %d/%d\n\n', studyName{1}, i_task, noTasks, i_study,length(studies)) 
            save([subFolderName '/_studies'], 'studies')
        end
        if task.useROM
            studies(i_study).loopParameters{end+1} = 'noVecs';
            studies(i_study).loopParametersArr{end+1} = noVecsArr;
            studies(i_study).loopParameters{end+1} = 'basisROM';
            studies(i_study).loopParametersArr{end+1} = basisROMcell;
            studies(i_study).tasks = tasks(:);
        end
        save([subFolderName '/_studies'], 'studies')
    end
end
fprintf('\n\nTotal time spent on study "%s": %12f seconds\n', studyName{1}, toc(t_start_study))  

%% Perform post plotting
close all
for i_study = 1:numel(studies)  
    study = studies(i_study);
    for i = 1:numel(study.postPlot)
        figure(i)
        printResultsToTextFiles(study,study.postPlot(i))
        if isa(study.postPlot(i).addCommands,'function_handle')
            figure(i)
            study.postPlot(i).addCommands(study,i_study,studies)
        end
    end
end

