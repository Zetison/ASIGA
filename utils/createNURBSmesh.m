function varCol = createNURBSmesh(varCol, model, M, degree)

if ~isfield(varCol{1}, 'meshFile')
    varCol{1}.meshFile = ['createNURBSmesh_' model];
end
eval(['varCol = ' varCol{1}.meshFile '(varCol, M, degree);'])

for i = 1:numel(varCol) % assume coreMethod to be the same in all domains
    if isfield(varCol{i},'nurbs')
        varCol{i} = repeatKnots(varCol{i},varCol{1}.coreMethod);
        varCol{i} = degenerateIGAtoFEM(varCol{i},varCol{1}.coreMethod);
    end
end
