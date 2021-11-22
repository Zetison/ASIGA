function task = createNURBSmesh(task)

if ~isfield(task.msh, 'meshFile')
    task.msh.meshFile = ['createNURBSmesh_' task.misc.model];
end
eval(['task = ' task.msh.meshFile '(task);'])

task = repeatKnots(task);
task = degenerateIGAtoFEM(task);
for i = 1:numel(task.varCol)
    if task.msh.explodeNURBS
        task.varCol{i}.nurbs = explodeNURBS(task.varCol{i}.nurbs);
    end
    if ~task.varCol{1}.boundaryMethod
        if i == 1
            task = defineDomains(task);
        else
            task.varCol = copySet(task.varCol,i-1, 'inner', 'innerCoupling');
            task.varCol = copySet(task.varCol,i, 'outer', 'outerCoupling');
        end
    end
end
PMLpatchFound = false;
for i = 1:numel(task.varCol{1}.nurbs)
    if isfield(task.varCol{1}.nurbs{i},'isPML') && any(task.varCol{1}.nurbs{i}.isPML)
        PMLpatchFound = true;
    end
end
if strcmp(task.misc.method,'PML') && ~PMLpatchFound
    task = createPML(task);
end
if ~task.varCol{1}.boundaryMethod
    for i = 1:numel(task.varCol)
        task.varCol = findCartesianAlignedBdry(task.varCol,i);
    end
end

function varCol = findCartesianAlignedBdry(varCol,domain)

connection = varCol{domain}.geometry.topology.connection;
names = {'yz','xz','xy'};
for j = 1:3
    patches = [];
    faces = [];
    for i = 1:numel(connection)
        at = zeros(2,3,'logical');
        patch = connection{i}.Attributes.master;
        midx = connection{i}.Attributes.midx;
        at(midx) = true;
        nurbs = subNURBS(varCol{domain}.nurbs(patch),'at',at.');
        Eps = 1e7*eps;
        if all(abs(nurbs{1}.coeffs(j,:)) < Eps)
            faces = [faces, midx];
            patches = [patches, patch];
        end
    end
    if ~isempty(faces)
        topset.Attributes.name = names{j};
        topset.Attributes.type = 'face';
        topset.Attributes.normal = 'outward';
        topset.item = cell(1,numel(patches));
        for i = 1:numel(patches)
            topset.item{i}.Text = faces(i);
            topset.item{i}.Attributes.patch = patches(i);
        end
        varCol{domain}.geometry.topologysets.set{end+1} = topset;
    end
end


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
    nurbsPML(i) = subNURBS(task.varCol{1}.nurbs(patch),'at',at.');
end
nurbsPML = normalBasedSurface2volume(nurbsPML,task.pml.t);
nurbsPML = makeUniformNURBSDegree(nurbsPML,task.varCol{1}.nurbs{1}.degree);
noNewKnots = task.varCol{1}.refinement(task.msh.M);
nurbsPML = insertKnotsInNURBS(nurbsPML,[0,0,noNewKnots(4)]);
for i = 1:numel(nurbsPML)
    nurbsPML{i}.isPML = [false,false,true];
end
task.varCol{1}.nurbs = uniteNURBS({task.varCol{1}.nurbs,nurbsPML});
task = repeatKnots(task);
task = degenerateIGAtoFEM(task);
task = addPMLtopology(task);


% task.varCol = copySet(task.varCol,1, 'Gamma_a', 'Gamma_a_PML');
% topset = task.varCol{1}.geometry.topologysets.set;
% idx = findSet(topset,'Gamma_a_PML');
% 
% connection = task.varCol{1}.geometry.topology.connection;
% item = task.varCol{1}.geometry.topologysets.set{idx}.item;
% for i = 1:numel(connection)
%     for j = 1:numel(item)
%         if connection{i}.Attributes.master == item{j}.Attributes.patch && connection{i}.Attributes.midx == item{j}.Text
%             item{j}.Attributes.patch = connection{i}.Attributes.slave;
%             item{j}.Text = connection{i}.Attributes.sidx;
%         end
%     end
% end
% task.varCol{1}.geometry.topologysets.set{idx}.item = item;






