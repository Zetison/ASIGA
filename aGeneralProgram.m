close all
clear all

addpath IGAfunctions
addpath NURBS
addpath SEM
addpath(genpath('NURBSgeometries'))
addpath otherFunctions
addpath(genpath('e3Dss'))
addpath(genpath('postProcessing'))
addpath integration
addpath subRoutines
addpath export_fig
% return

getTasks
runTasksInParallel = 0;
startMatlabPool
% Create folders
folderName = 'results';
if ~exist(folderName, 'dir')
    mkdir(folderName);
end

folderName = 'results/_studies';
if ~exist(folderName, 'dir')
    mkdir(folderName);
end
subFolderName = [folderName '/studies1'];
i = 1;
while exist(subFolderName, 'dir')
    i = i + 1;
    subFolderName = [folderName '/studies' num2str(i)];
end
mkdir(subFolderName);

for study_i = 1:numel(studies)    
    loopParameters = studies(study_i).loopParameters;
    loopParametersArr = studies(study_i).loopParametersArr;
    noTasks = length(studies(study_i).tasks);
    %% Allocate memory (in case of runTasksInParallel
    tasks = studies(study_i).tasks;
    if runTasksInParallel
%         for i_task = 1:length(tasks)
        parfor i_task = 1:length(tasks)
            tasks(i_task).task = aGeneralProgram_sub(tasks(i_task).task,loopParameters,runTasksInParallel);
%             [cond_number(i_task),p(i_task,:),abs_p(i_task,:),TS(i_task,:),p_ref(i_task,:),abs_p_ref(i_task,:),TS_ref(i_task,:),error_pAbs(i_task,:),error_p(i_task,:),energyError(i_task),L2Error(i_task),H1Error(i_task),surfaceError(i_task), ...
%                 k(i_task), f(i_task), dofs(i_task), nepw(i_task), tot_time(i_task), N(i_task), h_max(i_task), timeBuildSystem(i_task), analyticSolutionExist(i_task)] = aGeneralProgram_sub(tasks(i_task).task,loopParameters,runTasksInParallel);
            fprintf('\nCompleted task %d/%d in study %d/%d\n\n', i_task, noTasks, study_i,length(studies)) 
        end
        studies(study_i).tasks = tasks;
        studies(study_i).subFolderName = subFolderName;
        save([subFolderName '/_studies'], 'studies', '-v7.3')
    else
        for i_task = 1:noTasks
            task = tasks(i_task,1).task;
            aGeneralProgram_sub
            if plot3Dgeometry || plot2Dgeometry
                return
            end
            if useROM
                computeROMsolution
            end
            studies(study_i).tasks(i_task,1).task = task;
            studies(study_i).subFolderName = subFolderName;
            fprintf('\nCompleted task %d/%d in study %d/%d\n\n', i_task, noTasks, study_i,length(studies)) 
            save([subFolderName '/_studies'], 'studies', '-v7.3')
        end
        if useROM
            studies.loopParameters{end+1} = 'noVecs';
            studies.loopParametersArr{end+1} = noVecsArr;
            studies.loopParameters{end+1} = 'basisROM';
            studies.loopParametersArr{end+1} = basisROMcell;
            studies(study_i).tasks = tasks(:);
        end
        save([subFolderName '/_studies'], 'studies', '-v7.3')
    end
end

% plotResults

