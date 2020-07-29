variables = whos;
if ~exist('protectedVariables','var')
    protectedVariables = {};
end
for i = 1:length(variables)
    if ~any(ismember(variables(i).name,protectedVariables))
        container.(variables(i).name) = eval(variables(i).name);
    end
end