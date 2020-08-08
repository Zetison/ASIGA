function nurbs = getFreeNURBS(nurbs)

subnurbs = subNURBS(nurbs,'outwardPointingNormals', true);
subnurbs = explodeNURBS(subnurbs);
noSurfPatches = numel(subnurbs);
counter = 1;
for i = 1:noSurfPatches
    isFree = true;
    for j = 1:noSurfPatches
        if i ~= j && NURBSisEqual(subnurbs{i},subnurbs{j},1)
            isFree = false;
            break
        end
    end
    if isFree && ~NURBSisDegenerate(subnurbs{i})
        nurbs{counter} = subnurbs{i};
        counter = counter + 1;
    end
end
nurbs(counter:end) = [];