function geometry = getTopology(nurbs)

geometry = [];
if isempty(nurbs)
    return
end

Eps = 1e-10;
noPatches = numel(nurbs);
keys = cell(noPatches,1);
subnurbs = cell(noPatches,1);
for i = 1:noPatches
    subnurbs{i} = subNURBS(nurbs(i));
    keys{i} = zeros(1,numel(subnurbs{i}),'int32');
    for ii = 1:numel(subnurbs{i})
        coeffs = subnurbs{i}{ii}.coeffs;
        d_p = subnurbs{i}{ii}.d_p;
        d = subnurbs{i}{ii}.d;
        indices = 1:d_p;
        compScaling = [1,2,3]; % Add more uniquness to the key (due to symmetries in many geometries)
        compSum = sum(abs(subnurbs{i}{ii}.coeffs(1:d,:)),2);
        bdryMeasure = dot(compSum/norm(compSum),compScaling/norm(compScaling));
        keys{i}(ii) = int32(bdryMeasure*1e6);
        % Skip adding bdries that have zero measure
        for iii = 1:d_p
            zeroMeasure = true;
            temp = reshape(permute(coeffs,[1,iii+1,setdiff(indices,iii)+1]),d+1,subnurbs{i}{ii}.number(iii),[]);

            for l = 1:size(temp,3)
                temp2 = temp(1:d,:,l).';
                uniqueCoeffs = uniquetol(temp2,Eps,'ByRows',true, 'DataScale',max(norm2(temp2)), 'OutputAllIndices', true);
                if size(uniqueCoeffs,1) ~= 1
                    zeroMeasure = false;
                    break
                end
            end
            if zeroMeasure
                break
            end
        end
        if zeroMeasure
            keys{i}(ii) = -1;
        end
    end
end

counter = 1;
counterFree = 1;
freeBdries = zeros(noPatches,2);
for i = 1:noPatches
    if nurbs{i}.d_p == 0
        continue
    end
    for ii = 1:numel(subnurbs{i})
        if keys{i}(ii) < 0
            continue
        end
        % Add face with nonzero measure
        freeBdry = true;
        connectionFound = false;
        for j = 1:noPatches
            if ~ismember(keys{i}(ii),keys{j})
                continue
            end
            for jj = 1:numel(subnurbs{j})
                [matchingBdries,orient] = NURBSisEqual(subnurbs{i}{ii},subnurbs{j}{jj},1);
                if matchingBdries && ~(i == j && ii == jj)
                    freeBdry = false;
                end
                if ~(i == j && ii <= jj) && matchingBdries
                    if ~connectionFound
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
        end
        if freeBdry
            freeBdries(counterFree,:) = [i,ii];
            counterFree = counterFree + 1;
        end
    end
end
freeBdries(counterFree:end,:) = [];
if isempty(geometry)
    return
end
geometry.topologysets.set = {};
if ~isempty(freeBdries)
    % Find connected surfaces
    nurbsFaces = subNURBS(nurbs(freeBdries(:,1)),'at',num2cell(freeBdries(:,2).'));
    geometryBdry = getTopology(nurbsFaces);
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
    avgLen = zeros(1,numel(connMap));
    for i = 1:numel(connMap)
        nurbsPart = nurbsFaces(connMap{i});
        noCpts = 0;
        center = zeros(d,1);
        for j = 1:numel(nurbsPart)
            center = center + sum(nurbsPart{j}.coeffs(1:d,:),2);
            noCpts = noCpts + prod(nurbsPart{j}.number);
        end
        center = center/noCpts;
        noCpts = 0;
        for j = 1:numel(nurbsPart)
            avgLen(i) = avgLen(i) + sum(norm2((nurbsPart{j}.coeffs(1:d,:)-center).'));
            noCpts = noCpts + prod(nurbsPart{j}.number);
        end
        avgLen(i) = avgLen(i)/noCpts;
    end
    [~,I_outer] = max(avgLen);
    
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