%% Regression tests
% A successful run of this script should be performed prior to any push to git
startup
Eps = 1e-5;
testFolder = 'examples/taskFiles/tests/';
studyName = {'Venas2018iao_Figure6','Venas2020asi_Figure21','Venas2018iao_Table2','Venas2020asi_Figure7','Venas2019asi_Figure4','Venas2019asi_FigureB10B11B12',...
             'Venas2019asi_FigureA2A3'};
% studyName = {'Venas2018iao_Table2'};
stringShift = 60;
noFailedTests = 0;
for M_0 = 1:2 % Most test are only available for M_0 = 1 and M_0 = 2
    for i = 1:numel(studyName)
        fprintf(['\n%-' num2str(stringShift) 's'], ['Running test ''' studyName{i} ''' (M_0 = ' num2str(M_0) ') ...'])
        testFailed = false;
        try
            studiesCol = main(studyName{i},true,M_0);
            studiesCol_ref = load([testFolder studyName{i} '_M_0_' num2str(M_0) '.mat'],'studiesCol');
            studiesCol_ref = studiesCol_ref.studiesCol;
            for i_study = 1:numel(studiesCol{1})
                for i_task = 1:numel(studiesCol{1}(i_study).tasks)  
                    results = studiesCol{1}(i_study).tasks(i_task).task.results;
                    results_ref = studiesCol_ref{1}(i_study).tasks(i_task).task.results;
                    fieldNames = fieldnames(results_ref);
                    for field = fieldNames.'
                        entry = results.(field{1});
                        entry_ref = results_ref.(field{1});
                        if any(isnan(entry_ref(:)))
                            continue
                        end
                        reg_error = norm(entry(:)-entry_ref(:))/norm(entry_ref(:));
                        if reg_error > Eps
                            testFailed = true;
                            break
                        end
                    end
                    if testFailed
                        break
                    end
                end
                if testFailed
                    break
                end
            end
            if testFailed
                fprintf('test failed due to incorrect results (i_study=%d, i_task=%d)! (relative error = %g)', i_study, i_task, reg_error)
                noFailedTests = noFailedTests + 1;
            else
                fprintf('successfully!')
            end
        catch ME
            rethrow(ME)
            fprintf('Test failed due to runtime error!')
            noFailedTests = noFailedTests + 1;
        end
    end
end
fprintf(['\n\nNumber of failed tests: ' num2str(noFailedTests) '\n'])


return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create and store tests
% studyName = 'Venas2018iao_Figure6'; % getTask_Venas2018iao_Figure6        % Test CCBIE BEM (Simpson case)
% studyName = 'Venas2020asi_Figure21'; % getTask_Venas2020asi_Figure21      % Test MFS
% studyName = 'Venas2018iao_Table2'; % getTask_Venas2018iao_Table2        % Test FEM implementations and ASI implementation
% studyName = 'Venas2020asi_Figure7'; % getTask_Venas2020asi_Figure7      % Test SEM implementation
% studyName = 'Venas2019asi_Figure4'; % getTask_Venas2019asi_Figure4      % Test KDT implementation
% studyName = 'Venas2019asi_FigureB10B11B12'; % getTask_Venas2019asi_FigureB10B11B12      % Test RT implementation
studyName = 'Venas2019asi_FigureA2A3'; % getTask_Venas2019asi_FigureA2A3      % Test IENSG implementation

studiesColFull = main(studyName,true,M_0);
studiesCol = cell(1);
for i_study = 1:numel(studiesColFull{1})
    for i_task = 1:numel(studiesColFull{1}(i_study).tasks)  
        studiesCol{1}(i_study).tasks(i_task).task.results = studiesColFull{1}(i_study).tasks(i_task).task.results;
    end
end
save([testFolder studyName  '_M_0_' num2str(M_0) '.mat'],'studiesCol')
