function [j,setFound] = findSet(set,setName)
setFound = false;
for j = 1:numel(set)
    if strcmp(set{j}.Attributes.name,setName)
        setFound = true;
        break
    end
end
if ~setFound
    j = NaN;
end