function [nurbs, maxLengths] = autoRefineNURBS(nurbs,connection,h_max,dirs)

d_p = nurbs{1}.d_p;
if nargin < 4
    dirs = 1:d_p;
end

%% Create an easily searchable topology map to find connected patches needing the same knot insertions to ensure continuity
noPatches = numel(nurbs);
topologyMap = cell(1,noPatches);
for patch = 1:noPatches
    topologyMap{patch}(2*d_p) = struct();
    for i = 1:numel(connection)
        if connection{i}.Attributes.master == patch
            midx = connection{i}.Attributes.midx;
            topologyMap{patch}(midx).slave  = connection{i}.Attributes.slave;
            topologyMap{patch}(midx).sidx   = connection{i}.Attributes.sidx;
            topologyMap{patch}(midx).orient = connection{i}.Attributes.orient;
        end
    end
end

%% Find all edge lengths in nurbs in the directions dirs
edgeLen = edgeLengths(nurbs,dirs);

%% Search through the topology map to find path constrained to the same refinements
maxLengths = zeros(noPatches,d_p);
patchChecked = false(noPatches,d_p);
paths = {};
for patch = 1:noPatches
    if all(patchChecked(patch,:))
        continue
    end
    [maxLengths, patchChecked,paths_patch] = searchNURBS(nurbs,topologyMap,maxLengths, patchChecked,dirs,patch,{});
    if ~isempty(paths_patch{1})
        paths{end+1} = paths_patch;
    end
end

%% Go though all paths and find the maximal lengths along all elements
maxLengths = cell(1,numel(paths));
[~,~,orientMap] = getOrientPerms(d_p-1);
for i = 1:numel(paths)
    paths_patch = paths{i};
    maxLengths{i} = cell(1,numel(paths_patch));
    for j = 1:numel(paths_patch)
        path = paths_patch{j};
        midx = path(1).midx;
        if mod(midx,2) % midx is odd
            mdir = (midx+1)/2;
        else
            mdir = midx/2;
        end
        dir2 = circshift(1:d_p,1-mdir); % The other direction at which to insert knots
        dir2 = dir2(2:end);
        maxLengths{i}{j} = cell(1,numel(dir2));

        [maxLengths{i}{j}{:}] = edgeLen{path(1).master}{dir2};
        orient = zeros(numel(path),1);
        for ii = 2:numel(path)
            orient(ii) = path(ii-1).orient;

            midx = path(ii).midx;
            if mod(midx,2) % sidx is odd
                sdir = (midx+1)/2;
            else
                sdir = midx/2;
            end
            dir2_sidx = circshift(1:d_p,1-sdir); % The other direction at which to insert knots
            dir2_sidx = dir2_sidx(2:end);

            maxLengths_next = cell(1,numel(dir2_sidx));
            [maxLengths_next{:}] = edgeLen{path(ii).master}{dir2_sidx};

            % Take orient into account
            for jj = ii:-1:2
                flips = orientMap{orient(jj)+1,1};
                for jjj = 1:numel(flips)
                    flipIdx = flips(jjj);
                    maxLengths_next{flipIdx} = flip(maxLengths_next{flipIdx});
                end
                maxLengths_next = maxLengths_next(orientMap{orient(jj)+1,2});
            end

            % Update max lengths
            for jj = 1:numel(maxLengths{i}{j})
                maxLengths{i}{j}{jj} = max(maxLengths{i}{j}{jj},maxLengths_next{jj});
            end
        end
    end
end

%% Refine patches such that the elemnt sizes (based on the edge lengths) are roughly h_max
patchRefined = false(noPatches,d_p);
for i = 1:numel(paths)
    paths_patch = paths{i};
    for j = 1:numel(paths_patch)
        path = paths_patch{j};
        orient = zeros(numel(path),1);
        for ii = 1:numel(path)
            patch = path(ii).master;
            newKnots = cell(1,d_p);
            midx = path(ii).midx;
            if ii > 1
                orient(ii) = topologyMap{patch}(path(ii-1).sidx).orient;
            end
            if mod(midx,2) % midx is odd
                mdir = (midx+1)/2;
            else
                mdir = midx/2;
            end
            dir2 = circshift(1:d_p,1-mdir); % The other direction at which to insert knots
            dir2 = dir2(2:end);
            maxLengths_next = maxLengths{i}{j};

            % Take orient into account
            for jj = 2:ii
                flips = orientMap{orient(jj)+1,1};
                for jjj = 1:numel(flips)
                    flipIdx = flips(jjj);
                    maxLengths_next{flipIdx} = flip(maxLengths_next{flipIdx});
                end
                maxLengths_next = maxLengths_next(orientMap{orient(jj)+1,2});
            end
            for jj = 1:numel(dir2)
                if ~patchRefined(patch,dir2(jj))
                    newKnots{dir2(jj)} = insertNonUniform(nurbs{patch}.knots{dir2(jj)}, max(round(maxLengths_next{jj}/h_max)-1,0));
                end
            end
            nurbs(patch) = insertKnotsInNURBS(nurbs(patch),newKnots);
            patchRefined(patch,dir2) = true;
        end
    end
end

function [maxLengths, patchChecked, paths] = searchNURBS(nurbs,topologyMap,maxLengths,patchChecked,dirs,patch,paths)

for dir = dirs
    if patchChecked(patch,dir)
        continue
    else
        patchChecked(patch,dir) = true;
    end
    for i = [2*dir-1,2*dir]
        sub_paths = [];
        sub_paths.master = patch;
        slave = topologyMap{patch}(i).slave;
        sidx = topologyMap{patch}(i).sidx;
        if mod(i,2)
            sub_paths.midx = i;
            sub_paths.slave = slave;
            sub_paths.sidx = sidx;
            sub_paths.orient = topologyMap{patch}(i).orient;
        end
        if isempty(slave)
            if mod(i,2)
                paths{end+1} = sub_paths;
            end
            continue
        end
        [maxLengths, patchChecked, sub_paths] = searchNURBSdir(nurbs,topologyMap,maxLengths,patchChecked,slave,sub_paths,sidx,i);
        if mod(i,2)
            paths{end+1} = sub_paths;
        else
            if numel(sub_paths) > 1 % extended path has been found
                paths{end} = [fliplr(sub_paths(2:end)),paths{end}];
            end
        end
    end
end


function [maxLengths, patchChecked, paths] = searchNURBSdir(nurbs,topologyMap,maxLengths,patchChecked,patch,paths,midx,i)
if isempty(patch)
    return
end
% Continue to the other side of the patch
midx_prev = midx;
if mod(midx,2) % sidx is odd
    midx = midx + 1;
    dir = midx/2;
else
    dir = midx/2;
    midx = midx - 1;
end
if patchChecked(patch,dir)
    return
end

% Find next slave
slave = topologyMap{patch}(midx).slave;
sidx = topologyMap{patch}(midx).sidx;

patchChecked(patch,dir) = true;
if mod(i,2)
    paths(end+1).master = patch;
    paths(end).midx = midx;
    paths(end).slave = slave;
    paths(end).sidx = sidx;
    paths(end).orient = topologyMap{patch}(midx).orient;
else
    paths(end+1).master = patch;
    paths(end).midx = midx_prev;
    paths(end).slave = topologyMap{patch}(midx_prev).slave;
    paths(end).sidx = topologyMap{patch}(midx_prev).sidx;
    paths(end).orient = topologyMap{patch}(midx_prev).orient;
end
[maxLengths, patchChecked, paths] = searchNURBSdir(nurbs,topologyMap,maxLengths,patchChecked,slave,paths,sidx,i);

