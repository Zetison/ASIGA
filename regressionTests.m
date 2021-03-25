


M_0 = 1; % Most test are only availabel for M_0 = 1
Eps = eps;
testFolder = 'examples/taskFiles/tests/';
studyName = {'Venas2018iao_Figure6','Venas2020asi_Figure21','Venas2018iao_Table2'};
stringShift = 55;
noFailedTests = 0;
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
                for field = fieldNames
                    entry = results.(field{1});
                    entry_ref = results_ref.(field{1});
                    if any(isnan(entry_ref(:)))
                        continue
                    end
                    if norm(entry(:)-entry_ref(:))/norm(entry_ref) > Eps
                        testFailed = true;
                    end
                end
            end
        end
        if testFailed
            fprintf('test failed due to incorrect results!\n')
            noFailedTests = noFailedTests + 1;
        else
            fprintf('with success!')
        end
    catch
        fprintf('Test failed due to runtime error!\n')
        noFailedTests = noFailedTests + 1;
    end
end
fprintf(['\n\nNumber of failed tests: ' num2str(noFailedTests) '\n'])


return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create and store tests
% studyName = 'Venas2018iao_Figure6'; % getTask_Venas2018iao_Figure6
studyName = 'Venas2020asi_Figure21'; % getTask_Venas2020asi_Figure21
studyName = 'Venas2018iao_Table2'; % getTask_Venas2018iao_Table2

studiesColFull = main(studyName,true,M_0);
studiesCol = cell(1);
for i_study = 1:numel(studiesColFull{1})
    for i_task = 1:numel(studiesColFull{1}(i_study).tasks)  
        studiesCol{1}(i_study).tasks(i_task).task.results = studiesColFull{1}(i_study).tasks(i_task).task.results;
    end
end
save([testFolder studyName  '_M_0_' num2str(M_0) '.mat'],'studiesCol')
