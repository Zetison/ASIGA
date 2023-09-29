function task = createNURBSmesh(task)

if ~isfield(task.msh, 'meshFile')
    task.msh.meshFile = ['createNURBSmesh_' task.misc.model];
end
eval(['task = ' task.msh.meshFile '(task);'])

task = repeatKnots(task);
if task.msh.explodeNURBS
    for i = 1:numel(task.varCol)
        task.varCol{i}.nurbs = explodeNURBS(task.varCol{i}.nurbs);
    end
end
for i = 1:numel(task.varCol)
    if i == 1
        task = defineDomains(task);
    else
        task.varCol = copySet(task.varCol,i-1, 'inner', 'innerCoupling');
        task.varCol = copySet(task.varCol,i, 'outer', 'outerCoupling');
    end
end
if strcmp(task.misc.method,'PML') 
    PMLpatchFound = false;
    for i = 1:numel(task.varCol{1}.nurbs)
        if isfield(task.varCol{1}.nurbs{i},'isPML') && any(task.varCol{1}.nurbs{i}.isPML)
            PMLpatchFound = true;
        end
    end
    if ~PMLpatchFound
        if ~(isfield(task.pml,'refinement') || isfield(task.varCol{1},'refLength'))
            error('The refinement field must be specified using PML')
        end
        task = createPML(task);
    end
    task = addPMLtopology(task);
end

if task.msh.nonLinearParam
    task = convertAllToNonLinearParam(task);
end
if task.msh.autoRefine
    task.varCol = autoRefineDomains(task.varCol,task.msh.h_max,task.msh.dirs);
end
task = degenerateIGAtoFEM(task);

for i = 1:numel(task.varCol)
    task.varCol = findCartesianAlignedBdry(task.varCol,i);
end


function task = convertAllToNonLinearParam(task)

for i = 1:numel(task.varCol)
    task.varCol{i}.nurbs = convertToNonLinearParam(task.varCol{i}.nurbs);
end

function varCol = findCartesianAlignedBdry(varCol,domain)

connection = varCol{domain}.geometry.topology.connection;
names = {'yz','xz','xy','xy1'};
for j = 1:numel(names)
    patches = [];
    faces = [];
    for i = 1:numel(connection)
        patch = connection{i}.Attributes.master;
        midx = connection{i}.Attributes.midx;
        at = zeros(2,3,'logical');
        at(midx) = true;
        nurbs = subNURBS(varCol{domain}.nurbs(patch),'at',at.');
        Eps = 1e9*eps;
        if j == 4
            if all(abs(nurbs{1}.coeffs(3,:)) < Eps) && all(nurbs{1}.coeffs(1:3,:) > -Eps,'all')
                faces = [faces, midx];
                patches = [patches, patch];
            end
        else
            if all(abs(nurbs{1}.coeffs(j,:)) < Eps)
                faces = [faces, midx];
                patches = [patches, patch];
            end
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
nurbsPML = normalBasedSurface2volume(nurbsPML,task.pml.t,task.pml.X_bApprox);
nurbsPML = makeUniformNURBSDegree(nurbsPML,task.varCol{1}.nurbs{1}.degree);
if isfield(task.pml,'refinement')
    if nargin(task.pml.refinement) > 1
        noNewKnots = task.pml.refinement(task.msh.M,task.pml.t);
    else
        noNewKnots = task.pml.refinement(task.msh.M);
    end
    nurbsPML = insertKnotsInNURBS(nurbsPML,[0,0,noNewKnots]);
end
for i = 1:numel(nurbsPML)
    nurbsPML{i}.isPML = [false,false,true];
end
task.varCol{1}.nurbs = uniteNURBS({task.varCol{1}.nurbs,nurbsPML});
task = repeatKnots(task);
if task.msh.explodeNURBS
    for i = 1:numel(task.varCol)
        task.varCol{i}.nurbs = explodeNURBS(task.varCol{i}.nurbs);
    end
end


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






