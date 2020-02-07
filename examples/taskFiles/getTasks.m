function studies = getTasks(studyName)

studies = cell(0,1);

getDefaultTaskValues

counter = 1;
eval(['getTask_' studyName])
