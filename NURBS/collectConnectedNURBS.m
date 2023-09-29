function [nurbsCol,map] = collectConnectedNURBS(nurbs)

patchTop = getPatchTopology(nurbs);
noPatches = numel(nurbs);
indices = 1:noPatches;
nurbsCol = cell(1,noPatches);
map = cell(1,noPatches);
i = 1;
while ~isempty(indices)
    I = indices(1);
    j = 1;
    while isnan(I)
        I = indices(j);
        j = j + 1;
    end
    indices2 = getIndices(patchTop,I,I);
    indices = setdiff(indices, indices2);
    nurbsCol{i} = nurbs(indices2);
    map{i} = indices2;
    i = i + 1;
end
nurbsCol(i:end) = [];
map(i:end) = [];

function indices = getIndices(patchTop,indices,I)
if nargin < 3
    I = 1;
end
for i = 1:size(patchTop{I},1)
    j = patchTop{I}(i,1);
    if ~ismember(j,indices) && ~isnan(j)
        indices = [indices, j];
        indices = getIndices(patchTop,indices,j);
    end
end

