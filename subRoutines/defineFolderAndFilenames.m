

% if length(alpha_s) > 1
%     aspect = 'S';
% else
%     aspect = num2str(round(alpha_s*180/pi, 15, 'significant'));
% end
% if length(beta_s) > 1
%     elevation = 'S';
% else
%     elevation = num2str(round(beta_s*180/pi, 15, 'significant'));
% end
% if strcmp(scatteringCase, 'Sweep')
%     frequency = 'S';
% else
%     frequency = num2str(f/1000);
% end
% saveName = [model '_' BC '_' scatteringCase '_A'  aspect '_E' elevation '_F' frequency];

saveName = model;
for i = 1:length(loopParameters)
    temp2 = loopParameters{i};
    temp = task.(loopParameters{i});
    if isnumeric(temp)
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
% task.saveName = saveName;