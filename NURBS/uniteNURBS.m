function nurbs = uniteNURBS(nurbsPatches)
noPatches = numel(nurbsPatches);
noSubPatches = zeros(noPatches,1);

for i = 1:noPatches
    noSubPatches(i) = numel(nurbsPatches{i});
end

nurbs = cell(1,sum(noSubPatches));
for i = 1:noPatches
    nurbs(sum(noSubPatches(1:i-1))+1:sum(noSubPatches(1:i-1))+noSubPatches(i)) = nurbsPatches{i};
end