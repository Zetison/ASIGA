function nurbs = extractFreeSurface(varargin)
options = struct();

nurbs_vol = varargin{1};
if nargin > 1
    if numel(varargin) > 2
        newOptions = varargin(2:end);
    else
        newOptions = varargin{2};
    end
    options = updateOptions(options,newOptions);
end

if isfield(options,'topologysets')
    topologysets = options.topologysets;
else
    geometry = getTopology(nurbs_vol);
    topologysets = geometry.topologysets;
end

nurbs = cell(1,6*numel(nurbs_vol));
counter = 1;
for i = 1:numel(topologysets.set{1}.item)
    patch = topologysets.set{1}.item{i}.Attributes.patch;
    midx = topologysets.set{1}.item{i}.Text;
    at = zeros(2,3);
    at(midx) = 1;
    nurbs(counter) = subNURBS(nurbs_vol(patch),'at',at.','outwardPointingNormals',true);
    counter = counter + 1;
end
nurbs(counter:end) = [];