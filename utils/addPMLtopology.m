function task = addPMLtopology(task)
geometry = getTopology(task.varCol{1}.nurbs);
topset = geometry.topologysets.set;
idx = findSet(topset,'outer');
setOuter = geometry.topologysets.set{idx};
geometry.topologysets = task.varCol{1}.geometry.topologysets;
geometry.topologysets.set{idx} = setOuter;

task.varCol{1}.geometry.topologysets = geometry.topologysets;
task.varCol{1}.geometry.topology = geometry.topology;
if task.pml.dirichlet
    task.varCol = copySet(task.varCol,1, 'outer', 'homDirichlet');
end