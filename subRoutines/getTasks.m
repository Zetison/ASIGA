function studies = getTasks(studyName)

counter = 1;
studies = cell(0,1);

getDefaultTaskValues

eval(studyName)
