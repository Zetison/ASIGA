function geometry = getTopology(nurbs)

geometry.topologysets.set{1}.Attributes.name = 'inner';
geometry.topologysets.set{1}.Attributes.type = 'face';
geometry.topologysets.set{2}.Attributes.name = 'outer';
geometry.topologysets.set{2}.Attributes.type = 'face';

noPatches = numel(nurbs);
counter = 1;
counterFreeSurfInner = 1;
counterFreeSurfOuter = 1;
for i = 1:noPatches
    subnurbs_i = subNURBS(nurbs(i),'outwardPointingNormals', true);
    for ii = 1:numel(subnurbs_i)
        freeFace = true;
        for j = i:noPatches
            subnurbs_j = subNURBS(nurbs(j),'outwardPointingNormals', true);
            for jj = 1:numel(subnurbs_j)
                if ~(i == j && ii <= jj) && NURBSisEqual(subnurbs_i{ii},subnurbs_j{jj},1)
                    geometry.topology.connection{counter}.Attributes.master = i;
                    geometry.topology.connection{counter}.Attributes.midx = ii;
                    geometry.topology.connection{counter}.Attributes.slave = j;
                    geometry.topology.connection{counter}.Attributes.sidx = jj;
                    counter = counter + 1;
                    freeFace = false;
                end
            end
        end
        if freeFace
            if ii == 5
                geometry.topologysets.set{1}.item{counterFreeSurfInner}.Attributes.patch = i;
                geometry.topologysets.set{1}.item{counterFreeSurfInner}.Text = ii;
                counterFreeSurfInner = counterFreeSurfInner + 1;
            elseif ii == 6
                geometry.topologysets.set{2}.item{counterFreeSurfOuter}.Attributes.patch = i;
                geometry.topologysets.set{2}.item{counterFreeSurfOuter}.Text = ii;
                counterFreeSurfOuter = counterFreeSurfOuter + 1;
            end
        end
    end
end