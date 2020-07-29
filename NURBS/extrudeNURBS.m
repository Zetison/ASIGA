function nurbs = extrudeNURBS(varargin)
options = struct('extrudeDir', [0,0,0], 'flip', false);
nurbs = varargin{1};
if nargin > 1
    if numel(varargin) > 2
        newOptions = varargin(2:end);
    else
        newOptions = varargin{2};
    end
    options = updateOptions(options,newOptions);
end
nurbs = ensure3DNURBS(nurbs);
d = 3;
noPatches = numel(nurbs);
for patch = 1:noPatches
    coeffs = zeros([size(nurbs{patch}.coeffs), 2]);

    coeffs_i = nurbs{patch}.coeffs;
    subs = {[repmat({':'},1,ndims(coeffs)-1),{1}]};
    coeffs = subasgnArr(coeffs,coeffs_i,subs);
    subs = {[{1:d},repmat({':'},1,ndims(coeffs))]};
    coeffs_i = subasgnArr(coeffs_i,slc(coeffs_i,1:d,1)+options.extrudeDir.',subs);
    subs = {[repmat({':'},1,ndims(coeffs)-1),{2}]};
    coeffs = subasgnArr(coeffs,coeffs_i,subs);
    Zeta = [0,0,1,1];
    nurbs(patch) = createNURBSobject(coeffs,[nurbs{patch}.knots(:)', {Zeta}]);
end
if options.flip
    dir = nurbs{1}.d_p;
    nurbs = flipNURBSparametrization(nurbs,dir);
end
