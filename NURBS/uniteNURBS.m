function nurbs = uniteNURBS(nurbsPatches)
% Unites a set of sets of NURBS patches
noDomains = numel(nurbsPatches);
noSubPatches = zeros(noDomains,1);

for i = 1:noDomains
    noSubPatches(i) = numel(nurbsPatches{i});
end

nurbs = cell(1,sum(noSubPatches));
for i = 1:noDomains
    nurbs(sum(noSubPatches(1:i-1))+1:sum(noSubPatches(1:i-1))+noSubPatches(i)) = nurbsPatches{i};
end