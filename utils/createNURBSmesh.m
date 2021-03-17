function task = createNURBSmesh(task)

if ~isfield(task.msh, 'meshFile')
    task.msh.meshFile = ['createNURBSmesh_' task.misc.model];
end
eval(['task = ' task.msh.meshFile '(task);'])

task = repeatKnots(task);
task = degenerateIGAtoFEM(task);
for i = 1:numel(task.varCol) % assume coreMethod to be the same in all domains
    if ~task.varCol{1}.boundaryMethod
        if i == 1
            task.varCol = copySet(task.varCol,1, 'inner', 'Gamma');
            task.varCol = copySet(task.varCol,1, 'outer', 'Gamma_a');
            task.varCol = copySet(task.varCol,1, 'inner', 'Neumann');
            if numel(task.varCol) > 1
                task.varCol = copySet(task.varCol,2, 'outer', 'Gamma');
                task.varCol = copySet(task.varCol,2, 'outer', 'Neumann');
            end
        else
            task.varCol = copySet(task.varCol,i-1, 'inner', 'innerCoupling');
            task.varCol = copySet(task.varCol,i, 'outer', 'outerCoupling');
        end
    end
end
if strcmp(task.misc.method,'PML')
    task = createPML(task);
end

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


function task = createPML(task)
topset = task.varCol{1}.geometry.topologysets.set;

idx = findSet(topset,'Gamma_a');
item = topset{idx}.item;
nurbsPML = cell(1,numel(item));
for i = 1:numel(item)
    at = zeros(2,3,'logical');
    patch = topset{idx}.item{i}.Attributes.patch;
    midx = topset{idx}.item{i}.Text;
    at(midx) = true;
    nurbs = subNURBS(task.varCol{1}.nurbs(patch),'at',at.');
    nurbsPML(i) = normalBasedSurface2volume(nurbs{1},task.pml.t);
end
nurbsPML = makeUniformNURBSDegree(nurbsPML,task.varCol{1}.nurbs{1}.degree);
noNewKnots = task.varCol{1}.refinement(task.msh.M);
nurbsPML = insertKnotsInNURBS(nurbsPML,[0,0,noNewKnots(4)]);
for i = 1:numel(nurbsPML)
    nurbsPML{i}.isPML = true;
end
task.varCol{1}.nurbs = uniteNURBS({task.varCol{1}.nurbs,nurbsPML});
geometry = getTopology(task.varCol{1}.nurbs);
topset = geometry.topologysets.set;
idx = findSet(topset,'outer');
setOuter = geometry.topologysets.set{idx};
geometry.topologysets = task.varCol{1}.geometry.topologysets;
geometry.topologysets.set{idx} = setOuter;

task.varCol{1}.geometry = geometry;
if task.pml.dirichlet
    task.varCol = copySet(task.varCol,1, 'outer', 'homDirichlet');
end
task.varCol = copySet(task.varCol,1, 'Gamma_a', 'Gamma_a_PML');
topset = task.varCol{1}.geometry.topologysets.set;
idx = findSet(topset,'Gamma_a_PML');

connection = task.varCol{1}.geometry.topology.connection;
item = task.varCol{1}.geometry.topologysets.set{idx}.item;
for i = 1:numel(connection)
    for j = 1:numel(item)
        if connection{i}.Attributes.master == item{j}.Attributes.patch && connection{i}.Attributes.midx == item{j}.Text
            item{j}.Attributes.patch = connection{i}.Attributes.slave;
            item{j}.Text = connection{i}.Attributes.sidx;
        end
    end
end
task.varCol{1}.geometry.topologysets.set{idx}.item = item;






