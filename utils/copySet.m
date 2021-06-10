function varCol = copySet(varCol,domain, source, target)
if ~isfield(varCol{domain},'geometry')
    varCol{domain}.geometry = getTopology(varCol{domain}.nurbs);
end
topset = varCol{domain}.geometry.topologysets.set;
idx = findSet(topset,source);
[~,setFound] = findSet(topset,target);
if ~setFound
    topset{end+1} = topset{idx};
    topset{end}.Attributes.name = target;
end
varCol{domain}.geometry.topologysets.set = topset;
