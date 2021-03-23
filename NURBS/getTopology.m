function geometry = getTopology(nurbs)

noPatches = numel(nurbs);
counter = 1;
counterFreeSurfInner = 1;
counterFreeSurfOuter = 1;
inneritem = cell(0);
outeritem = cell(0);
for i = 1:noPatches
    subnurbs_i = subNURBS(nurbs(i));
    for ii = 1:numel(subnurbs_i)
        % Skip adding faces that has zero measure
        coeffs = subnurbs_i{ii}.coeffs(1:3,:,:);
        [~,n,m] = size(coeffs);
        zeroMeasure = true;
        for l = 1:m
            uniqueCoeffs = uniquetol(coeffs(:,:,l).',1e7*eps,'ByRows',true, 'DataScale',max(norm2(coeffs(:,:).')), 'OutputAllIndices', true);
            if size(uniqueCoeffs,1) ~= 1
                zeroMeasure = false;
                break
            end
        end
        if zeroMeasure
            continue
        else
            zeroMeasure = true;
            coeffs = permute(coeffs,[1,3,2]);
            for l = 1:n
                uniqueCoeffs = uniquetol(coeffs(:,:,l).',1e7*eps,'ByRows',true, 'DataScale',max(norm2(coeffs(:,:).')), 'OutputAllIndices', true);
                if size(uniqueCoeffs,1) ~= 1
                    zeroMeasure = false;
                    break
                end
            end
        end
        if zeroMeasure
            continue
        end
        
        % Add face with nonzero measure
        freeFace = true;
        connectionFound = false;
        for j = 1:noPatches
            subnurbs_j = subNURBS(nurbs(j));
            for jj = 1:numel(subnurbs_j)
                matchingFaces = NURBSisEqual(subnurbs_i{ii},subnurbs_j{jj},1);
                if matchingFaces && ~(i == j && ii == jj)
                    freeFace = false;
                end
                if ~(i == j && ii <= jj) && matchingFaces
                    if ~connectionFound
                        geometry.topology.connection{counter}.Attributes.master = i;
                        geometry.topology.connection{counter}.Attributes.midx = ii;
                        geometry.topology.connection{counter}.Attributes.slave = j;
                        geometry.topology.connection{counter}.Attributes.sidx = jj;
                        connectionFound = true;
                        counter = counter + 1;
                    end
                end
            end
        end
        if freeFace
            if ii == 6
                outeritem{counterFreeSurfOuter}.Attributes.patch = i;
                outeritem{counterFreeSurfOuter}.Text = ii;
                counterFreeSurfOuter = counterFreeSurfOuter + 1;
            else
                inneritem{counterFreeSurfInner}.Attributes.patch = i;
                inneritem{counterFreeSurfInner}.Text = ii;
                counterFreeSurfInner = counterFreeSurfInner + 1;
            end
        end
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


