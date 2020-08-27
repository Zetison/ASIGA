function studiesCol = getTasks(studyName)

studiesCol = cell(1,numel(studyName));

for i_coll = 1:numel(studyName)
    eval(['studiesCol{i_coll} = getTask_' studyName{i_coll}])
end

