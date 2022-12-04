function studiesCol = main(studyName,runRegressionTests,M_0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Main function of ASIGA
% Author: Jon Vegard Venås
% E-mail: JonVegard.Venas@sintef.no
% Institute: SINTEF Digital
% Release: 3
% Release date: 23/07/2022

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
            progressBars = studiesCol{i_col}(i_study).tasks(1).task.misc.progressBars;
            nProgressStepSize = ceil(noTasks/1000);
            if progressBars
                ppm = ParforProgMon('Simulating tasks: ', noTasks);
            else
                ppm = NaN;
            end
            tasks = studiesCol{i_col}(i_study).tasks(i_task);
            parfor i_task = 1:noTasks
                if progressBars && mod(i_task,nProgressStepSize) == 0
                    ppm.increment();
                end
                
                %% Run main subroutine
                tasks(i_task).task = main_sub(tasks(i_task).task,loopParameters,printLog,resultsFolder);
                if ~runRegressionTests
                    fprintf('\nCase %s: Completed task %d/%d in study %d/%d in %f seconds.\n\n', studyName{i_col}, i_task, noTasks, i_study, noStudies, tasks(i_task).task.varCol{1}.tot_time) 
                end
            end
            studiesCol{i_col}(i_study).tasks = tasks;
        else
            for i_task = 1:noTasks
                plot3Dgeometry      = studiesCol{i_col}(i_study).tasks(i_task).task.prePlot.plot3Dgeometry;
                abortAfterPlotting  = studiesCol{i_col}(i_study).tasks(i_task).task.prePlot.abortAfterPlotting;
                preProcessOnly      = studiesCol{i_col}(i_study).tasks(i_task).task.misc.preProcessOnly;

                %% Run main subroutine
                studiesCol{i_col}(i_study).tasks(i_task).task = main_sub(studiesCol{i_col}(i_study).tasks(i_task).task,loopParameters,printLog,resultsFolder);
                if (plot3Dgeometry || preProcessOnly) && abortAfterPlotting
                    return
                end
                if printLog
                    fprintf('\nCase %s: Completed task %d/%d in study %d/%d\n\n', studyName{i_col}, i_task, noTasks, i_study, noStudies) 
                end
            end
        end 
        studiesCol{i_col}(i_study).resultsFolder = resultsFolder;
        if studiesCol{i_col}(i_study).saveStudies
            save([resultsFolder '/studiesCol'], 'studiesCol')
        end
        study = studiesCol{i_col}(i_study);
        for i = 1:numel(study.postPlot)
            study.postPlot(i).plotResults = false;
            if study.postPlot(i).printResults
                printResultsToTextFiles(study,study.postPlot(i))
            end
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
                    xname = strrep(xname, 'varCol', 'd');
                    
                    xname = regexprep(xname,'[^a-zA-Z1-9\s]','');
                    yname = regexprep(yname,'[^a-zA-Z1-9\s]','');
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

