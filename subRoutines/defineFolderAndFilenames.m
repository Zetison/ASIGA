function saveName = defineFolderAndFilenames(task,loopParameters)
saveName = task.misc.model;
for i = 1:length(loopParameters)
    temp2 = loopParameters{i}(find(loopParameters{i} == '.',1,'last')+1:end);
    eval(['temp = task.' loopParameters{i} ';'])
    if ~ischar(temp)
        temp = num2str(temp);
    end
    if isstruct(temp)
        fieldNames = fieldnames(temp);
        temp2 = fieldNames{1};
        temp = temp.(temp2);
    end
    temp = erase(temp,':');
    temp = erase(temp,' ');
    switch loopParameters{i}
        case {'formulation','IEbasis','method','coreMethod','BC'}
            saveName = [saveName '_' temp];
        otherwise
            saveName = [saveName '_' temp2 temp];
    end
end