
saveName = model;
for i = 1:length(loopParameters)
    temp2 = loopParameters{i};
    temp = task.(loopParameters{i});
    if ~ischar(temp)
        temp = num2str(temp);
    end
    if isstruct(temp)
        fieldNames = fieldnames(temp);
        temp2 = fieldNames{1};
        temp = temp.(temp2);
    end
    switch loopParameters{i}
        case {'formulation','IEbasis','method','coreMethod','BC'}
            saveName = [saveName '_' temp];
        otherwise
            saveName = [saveName '_' temp2 temp];
    end
end