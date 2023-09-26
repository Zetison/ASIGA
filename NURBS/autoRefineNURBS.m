function [nurbs, maxLengths,patchTrulyRefined] = autoRefineNURBS(nurbs,connection,h_max,dirs,start_patch,newKnots)
% Assuming all patches are compatible, this algorithm automatically refines
% all patches in nurbs such that maximal element side length is roughly
% h_max
d_p = nurbs{1}.d_p; % Assume all patches have the same number of parametric directions
if nargin < 4
    dirs = 1:d_p;
end

%% Create an easily searchable topology map to find connected patches needing the same knot insertions to ensure continuity
noPatches = numel(nurbs);
if nargin < 5
    start_patch = 1;
end
if nargin < 6
    newKnots = {};
end
topologyMap = createTopologyMap(connection,noPatches,d_p);

%% Find all edge lengths in nurbs in the directions dirs
edgeLen = edgeLengths(nurbs,1:d_p);

%% Search through the topology map to find path constrained to the same refinements and do refinement
patchRefined = false(noPatches,d_p);
patchTrulyRefined = false(noPatches,d_p);
[~,~,orientMap] = getOrientPerms(d_p-1);
for patch = [start_patch,setdiff(1:noPatches,start_patch)]
    if all(patchRefined(patch,:))
        continue
    end
    for i = 1:numel(dirs)
        refDir = dirs(i);
        if patch == start_patch
            newKnotsStart = newKnots{i};
        else
            newKnotsStart = [];
        end
        [nurbs, patchRefined, ~, maxLengths, patchTrulyRefined] = searchNURBS(nurbs,topologyMap,orientMap,edgeLen,patchRefined,refDir,false,h_max,patch,newKnotsStart,patchTrulyRefined);
    end
end

function [nurbs, patchRefined, patchChecked, maxLengths, patchTrulyRefined] = searchNURBS(nurbs,topologyMap,orientMap,edgeLen,patchRefined,refDir,...
                                                                        refine,h_max,rootPatch,newKnotsStart,patchTrulyRefined,master,patchChecked,maxLengths)
% This function goes through all patches that must be refined if master is
% refined in direction refDir: The maximum element edges are first computed
% and then the patches are refined based on this such that the element
% edges in refDir has roughly length h_max

if nargin < 12
    master = rootPatch;
    patchChecked = patchRefined;
    noElems = numel(unique(nurbs{master}.knots{refDir}))-1;
    maxLengths = -Inf(1,noElems);
end
if patchChecked(master,refDir)
    return
end
d_p = nurbs{master}.d_p;
pathDirs = setdiff(1:d_p,refDir);
mIndices = [2*pathDirs-1,2*pathDirs];

patchChecked(master,refDir) = true;

if refine
    newKnots = cell(1,d_p);
    if isempty(newKnotsStart)
        newKnots{refDir} = insertNonUniform(nurbs{master}.knots{refDir}, max(round(maxLengths/h_max)-1,0));
    else
        newKnots{refDir} = newKnotsStart;
    end
    patchTrulyRefined(master,refDir) = numel(newKnots{refDir}) > 0;
    patchRefined(master,refDir) = true;
    nurbs(master) = insertKnotsInNURBS(nurbs(master),newKnots);
else
    maxLengths = max(maxLengths,edgeLen{master}{refDir});
end
for midx = mIndices % Loop through all branches in the tree from the node "patch"
    slave = topologyMap{master}(midx).slave;
    if isempty(slave) % master has an open boundary at edge midx
        continue
    end

    % Directions in master for the the master-slave interface
    mdir = ceil(midx/2);
    dir2_midx = circshift(1:d_p,1-mdir);
    dir2_midx = dir2_midx(2:end);

    % Directions in slave for the master-slave interface
    sidx = topologyMap{master}(midx).sidx;
    if mod(sidx,2) % sidx is odd
        sdir = (sidx+1)/2;
    else
        sdir = sidx/2;
    end
    dir2_sidx = circshift(1:d_p,1-sdir);
    dir2_sidx = dir2_sidx(2:end);

    % Find the matching direction in the slave patch to refine matching the master patch refinement direction
    orient = topologyMap{master}(midx).orient;
    indices = orientMap{orient+1,2};
    refDirSlave = dir2_sidx(indices(dir2_midx == refDir)); 

    % Take orient into account and flip knot vectors if necessary
    flips = orientMap{orient+1,1};
    flipKnotVector = ismember(refDirSlave,dir2_sidx(flips));
    if flipKnotVector 
        maxLengths = flip(maxLengths);
    end

    % Continue search/refinement in the slave patch
    [nurbs, patchRefined, patchChecked,maxLengths, patchTrulyRefined] ...
        = searchNURBS(nurbs,topologyMap,orientMap,edgeLen,patchRefined,refDirSlave,refine,h_max,rootPatch,newKnotsStart,patchTrulyRefined,...
                        slave,patchChecked,maxLengths);

    % Flip back the refinement direction
    if flipKnotVector
        maxLengths = flip(maxLengths);
    end
end
if master == rootPatch && ~refine % Max lengths have been computed and we are now ready to refine the patches in the tree
    refine = true;
    patchChecked = patchRefined; % reset patchChecked to go through the tree the same way again
    [nurbs, patchRefined, patchChecked,maxLengths, patchTrulyRefined] ...
        = searchNURBS(nurbs,topologyMap,orientMap,edgeLen,patchRefined,refDir,refine,h_max,rootPatch,newKnotsStart,patchTrulyRefined,...
                        master,patchChecked,maxLengths);
end
