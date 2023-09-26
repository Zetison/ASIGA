function [nurbs,surf2volMap] = extractFreeSurface(varargin)
% Return the NURBS for the first free surface of the volumetric patches
% (nurbs_vol) found in geometry.topologysets.set constructed from nurbs_vol
options = struct('topoSetName','');

nurbs_vol = varargin{1};
if nargin > 1
    if numel(varargin) > 2
        newOptions = varargin(2:end);
    else
        newOptions = varargin{2};
    end
    options = updateOptions(options,newOptions);
end

if numel(nurbs_vol) == 1
    nurbs = subNURBS(nurbs_vol,'outwardPointingNormals',true);
else
    if isfield(options,'topologysets')
        topologysets = options.topologysets;
    else
        geometry = getTopology(nurbs_vol);
        topologysets = geometry.topologysets;
    end
    
    if isempty(options.topoSetName)
        topset = topologysets.set{1};
    else
        idx = findSet(geometry.topologysets.set,options.topoSetName);
        topset = topologysets.set{idx};
    end
    nurbs = cell(1,6*numel(nurbs_vol));
    surf2volMap = zeros(numel(nurbs),2);
    counter = 1;
    for i = 1:numel(topset.item)
        patch = topset.item{i}.Attributes.patch;
        midx = topset.item{i}.Text;
        at = zeros(2,3);
        at(midx) = 1;
        nurbs(counter) = subNURBS(nurbs_vol(patch),'at',at.','outwardPointingNormals',true);
        counter = counter + 1;
        surf2volMap(counter,:) = [patch,midx];
    end
    nurbs(counter:end) = [];
    surf2volMap(counter:end,:) = [];
end