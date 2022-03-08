function studiesCol = main(studyName,runRegressionTests,M_0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Main function of ASIGA
% Author: Jon Vegard Venås
% E-mail: JonVegard.Venas@sintef.no
% Institute: SINTEF Digital
% Release: 2
% Release date: 29/07/2020

startup
if nargin < 1
    studyName = availableStudies();
end

if nargin < 2
    runRegressionTests = false;
    M_0 = NaN;
end

%% Extract tasks
if ~iscell(studyName)
    studyName = {studyName};
end
studiesCol = getTasks(studyName,runRegressionTests,M_0);


%% Perform studies
for i_col = 1:numel(studiesCol)
    t_start_study = tic;
    noStudies = numel(studiesCol{i_col}) ;
    for i_study = 1:noStudies    
        loopParameters = studiesCol{i_col}(i_study).loopParameters;
        loopParametersArr = studiesCol{i_col}(i_study).loopParametersArr;
        noTasks = length(studiesCol{i_col}(i_study).tasks);

        tasks = studiesCol{i_col}(i_study).tasks;
        resultsFolder = [folderName '/' studiesCol{i_col}(i_study).subFolderName];
        if ~exist(resultsFolder, 'dir')
            mkdir(resultsFolder)
        end
        diary([resultsFolder '/ASIGA.log'])

        noCoresToUse = studiesCol{i_col}(i_study).noCoresToUse;
        runTasksInParallel = studiesCol{i_col}(i_study).runTasksInParallel;
        startMatlabPool
        printLog = ~runTasksInParallel && ~runRegressionTests;
        if runTasksInParallel
            progressBars = tasks(1).task.misc.progressBars;
            nProgressStepSize = ceil(length(tasks)/1000);
            if progressBars
                ppm = ParforProgMon('Simulating tasks: ', length(tasks));
            else
                ppm = NaN;
            end
            parfor i_task = 1:length(tasks)
                if progressBars && mod(i_task,nProgressStepSize) == 0
                    ppm.increment();
                end
                tasks(i_task).task = main_sub(tasks(i_task).task,loopParameters,printLog,resultsFolder);
                if ~runRegressionTests
                    fprintf('\nCase %s: Completed task %d/%d in study %d/%d in %f seconds.\n\n', studyName{i_col}, i_task, noTasks, i_study, noStudies, tasks(i_task).task.varCol{1}.tot_time) 
                end
            end
            studiesCol{i_col}(i_study).tasks = tasks;
            studiesCol{i_col}(i_study).resultsFolder = resultsFolder;
        else
            for i_task = 1:noTasks
                tasks(i_task).task = main_sub(tasks(i_task).task,loopParameters,printLog,resultsFolder);
                if (tasks(i_task).task.prePlot.plot3Dgeometry || tasks(i_task).task.prePlot.plot2Dgeometry) && tasks(i_task).task.prePlot.abortAfterPlotting
                    return
                end
                if tasks(i_task).task.rom.useROM
                    basisROMcell = studiesCol{i_col}(i_study).basisROMcell;
                    omega_ROM = studiesCol{i_col}(i_study).omega_ROM;
                    noVecsArr = studiesCol{i_col}(i_study).noVecsArr;
                    tasks = computeROMsolution(tasks,i_task,basisROMcell,omega_ROM,noVecsArr,printLog);
                end
                if printLog
                    fprintf('\nCase %s: Completed task %d/%d in study %d/%d\n\n', studyName{i_col}, i_task, noTasks, i_study, noStudies) 
                end
            end
            studiesCol{i_col}(i_study).tasks = tasks;
            if tasks(i_task).task.rom.useROM
                studiesCol{i_col}(i_study).loopParameters{end+1} = 'rom.noVecs';
                studiesCol{i_col}(i_study).loopParametersArr{end+1} = noVecsArr;
                studiesCol{i_col}(i_study).loopParameters{end+1} = 'rom.basisROM';
                studiesCol{i_col}(i_study).loopParametersArr{end+1} = basisROMcell;
                studiesCol{i_col}(i_study).tasks = tasks(:);
            end
            studiesCol{i_col}(i_study).resultsFolder = resultsFolder;
        end 
        if studiesCol{i_col}(i_study).saveStudies
            save([resultsFolder '/studiesCol'], 'studiesCol')
        end
    end
    close all
    for i_study = 1:numel(studiesCol{i_col})  
        study = studiesCol{i_col}(i_study);
        for i = 1:numel(study.postPlot)
            if study.postPlot(i).plotResults || study.postPlot(i).printResults
                figure(i)
                printResultsToTextFiles(study,study.postPlot(i))
                if isa(study.postPlot(i).addCommands,'function_handle') && i_study == numel(studiesCol{i_col}) 
                    figure(i)
                    study.postPlot(i).addCommands(study,i_study,studiesCol{i_col})
                    if isempty(study.postPlot(i).subFolderName)
                        subFolderName = study.resultsFolder;
                    else
                        subFolderName = study.postPlot(i).subFolderName;
                    end
                    xname = study.postPlot(i).xname;
                    yname = study.postPlot(i).yname;
                    xname(xname == '.') = [];
                    yname(yname == '.') = [];
                    model = study.tasks(1).task.misc.model;
                    savefig([subFolderName '/plot_' model '_' yname 'VS' xname])
                end
            end
        end
    end 
    if printLog
        fprintf('\n\nTotal time spent on case "%s": %12f seconds\n', studyName{i_col}, toc(t_start_study)) 
    end
end

