function studies = getTasks(studyName)

studies = cell(0,1);

counter = 1;
for i = 1:numel(studyName)
    getDefaultTaskValues
    eval(['getTask_' studyName{i}])
end

