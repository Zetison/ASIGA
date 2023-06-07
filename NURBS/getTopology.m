function geometry = getTopology(nurbs)

geometry = [];
if isempty(nurbs)
    return
end

noPatches = numel(nurbs);

%% Extract all interfaces (subnurbs) from the patches and assign each of these a unique key (stored in keys)
keys = cell(noPatches,1);
subnurbs = cell(noPatches,1);
for i = 1:noPatches
    if nurbs{i}.d_p == 0
        continue
    end
    subnurbs{i} = subNURBS(nurbs(i));
    keys{i} = zeros(1,numel(subnurbs{i}),'int32');
    for ii = 1:numel(subnurbs{i})
        d_p = subnurbs{i}{ii}.d_p;
        d = subnurbs{i}{ii}.d;
        compScaling = 1:d; % Add more uniquness to the key (due to symmetries in many geometries)
        compSum = sum(subnurbs{i}{ii}.coeffs(1:d,:),2);
        bdryMeasure = dot(compSum,compScaling);
        keys{i}(ii) = int32(bdryMeasure/10^ceil(log10(abs(bdryMeasure)))*1e9);
        if d_p >= 1
            % Skip adding bdries that have zero measure
            if NURBShasZeroMeasure(subnurbs{i}{ii})
                keys{i}(ii) = intmax;
            end
        end
    end
end

%% Get topology and find free boundaries
counter = 1;
counterFree = 1;
freeBdries = zeros(noPatches,2);
for i = 1:noPatches
    if nurbs{i}.d_p == 0
        continue
    end
    for ii = 1:numel(subnurbs{i})
        if keys{i}(ii) == intmax
            continue
        end
        % Add face with nonzero measure
        freeBdry = true;
        connectionFound = false;
        for j = 1:noPatches
            for jj = find(keys{j} == keys{i}(ii))
                [matchingBdries,orient] = NURBSisEqual(subnurbs{i}{ii},subnurbs{j}{jj},1);
                if matchingBdries && ~(i == j && ii == jj)
                    freeBdry = false;
                end
                if ~(i == j && ii <= jj) && matchingBdries && ~connectionFound
                    geometry.topology.connection{counter}.Attributes.master = i;
                    geometry.topology.connection{counter}.Attributes.midx = ii;
                    geometry.topology.connection{counter}.Attributes.slave = j;
                    geometry.topology.connection{counter}.Attributes.sidx = jj;
                    geometry.topology.connection{counter}.Attributes.orient = orient;
                    connectionFound = true;
                    counter = counter + 1;
                end
            end
        end
        if freeBdry
            freeBdries(counterFree,:) = [i,ii];
            counterFree = counterFree + 1;
        end
    end
end

%% Create topology sets
freeBdries(counterFree:end,:) = [];
if isempty(geometry)
    return
end
geometry.topologysets.set = {};
if ~isempty(freeBdries)
    % Find connected surfaces
    nurbsFaces = subNURBS(nurbs(freeBdries(:,1)),'at',num2cell(freeBdries(:,2).'));
    geometryBdry = getTopology(nurbsFaces);
    if isempty(geometryBdry)
        return
    end
    connections = geometryBdry.topology.connection;
    noFaces = numel(nurbsFaces);
    connMap = num2cell(1:noFaces);
    for i = 1:numel(connections)
        master = connections{i}.Attributes.master;
        slave = connections{i}.Attributes.slave;
        for master_i = 1:noFaces
            if ismember(master,connMap{master_i})
                break
            end
        end
        for slave_i = 1:noFaces
            if ismember(slave,connMap{slave_i})
                break
            end
        end
        if master_i ~= slave_i
            connMap{master_i} = sort([connMap{master_i}, connMap{slave_i}]);
            connMap{slave_i} = [];
        end
    end
    
    % Find the outer surface
    connMap = connMap(~cellfun('isempty',connMap));
    d = nurbsFaces{1}.d;
    V_max = -Inf;
    for i = 1:numel(connMap)
        nurbsPart = nurbsFaces(connMap{i});
%         noCpts = 0;
        noCpts = zeros(numel(nurbsPart),1);
        for j = 1:numel(nurbsPart)
            noCpts(j) = prod(nurbsPart{j}.number);
        end

        % Collect all control points into a single array X
        X = zeros(sum(noCpts),d);
        counter = 1;
        for j = 1:numel(nurbsPart)
            X(counter:counter+noCpts(j)-1,:) = nurbsPart{j}.coeffs(1:d,:).';
            counter = counter + noCpts(j);
        end

        % Find the surface having the largest volume and define this to be
        % the outermost surface
        try
            [~, V] = convhull(X);
            if V > V_max
                I_outer = i;
                V_max = V;
            end
        catch
            warning('Could not compute the convex hull. The surface may not be closed.')
        end
    end
    
    % Generate topologysets
    counterFreeSurfOuter = 1;
    counterFreeSurfInner = 1;
    inneritem = cell(0);
    outeritem = cell(0);
    for i = 1:numel(nurbsFaces)
        if ismember(i,connMap{I_outer})
            outeritem{counterFreeSurfOuter}.Attributes.patch = freeBdries(i,1);
            outeritem{counterFreeSurfOuter}.Text = freeBdries(i,2);
            counterFreeSurfOuter = counterFreeSurfOuter + 1;
        else
            inneritem{counterFreeSurfInner}.Attributes.patch = freeBdries(i,1);
            inneritem{counterFreeSurfInner}.Text = freeBdries(i,2);
            counterFreeSurfInner = counterFreeSurfInner + 1;
        end
    end
    idx = 1;
    if ~isempty(inneritem)
        geometry.topologysets.set{idx}.item = inneritem;
        geometry.topologysets.set{idx}.Attributes.name = 'inner';
        geometry.topologysets.set{idx}.Attributes.type = 'face';
        geometry.topologysets.set{idx}.Attributes.normal = 'inward';
        idx = idx + 1;
    end
    if ~isempty(outeritem)
        geometry.topologysets.set{idx}.item = outeritem;
        geometry.topologysets.set{idx}.Attributes.name = 'outer';
        geometry.topologysets.set{idx}.Attributes.type = 'face';
        geometry.topologysets.set{idx}.Attributes.normal = 'outward';
    end
end