function varCol = autoRefineDomains(varCol,h_max,dirs)
noDomains = numel(varCol);
indices = cell(1,noDomains);
nurbs = [];
counter = 0;
d_p_max = -Inf;
for i = 1:noDomains
    noPatches = numel(varCol{i}.nurbs);
    indices{i} = counter+(1:noPatches);
    nurbs = [nurbs,varCol{i}.nurbs];
    counter = counter + noPatches;
    if varCol{i}.nurbs{1}.d_p > d_p_max
        d_p_max = varCol{i}.nurbs{1}.d_p;
    end
end
handleBoundaryMethodCase = false;
if noDomains > 1 && varCol{1}.nurbs{1}.d_p < d_p_max && noDomains < 4
    handleBoundaryMethodCase = true;
    nurbs = varCol{2}.nurbs;
end

geometry = getTopology(nurbs);
nurbs = autoRefineNURBS(nurbs,'connection',geometry.topology.connection,'h_max',h_max,'dirs',dirs(dirs <= d_p_max));
if handleBoundaryMethodCase
    if noDomains == 3
        varCol{3}.nurbs = subNURBS(nurbs,'at',[0,0;0,0;1,0]);
    end
    varCol{1}.nurbs = subNURBS(nurbs,'at',[0,0;0,0;0,1]);
else
    for i = 1:noDomains
        varCol{i}.nurbs = nurbs(indices{i});
    end
end