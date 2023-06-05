function topologyMap = createTopologyMap(connection,noPatches,d_p)
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