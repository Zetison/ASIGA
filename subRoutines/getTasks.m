function studiesCol = getTasks(studyName,runRegressionTests,M_0)

studiesCol = cell(1,numel(studyName));

for i_coll = 1:numel(studyName)
    if runRegressionTests
        eval(['studiesCol{i_coll} = getTask_' studyName{i_coll} '(M_0);'])
    else
        eval(['studiesCol{i_coll} = getTask_' studyName{i_coll} '();'])
    end
end

