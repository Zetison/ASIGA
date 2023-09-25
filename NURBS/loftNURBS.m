function nurbs = loftNURBS(nurbsCol,p,dir)
if nargin < 2
    p = 1;
end
if nargin < 3
    dir = 3;
end
for i = 1:numel(nurbsCol)
    nurbsCol{i} = ensure3DNURBS(nurbsCol{i});
end
noPatches = numel(nurbsCol{1});
n = numel(nurbsCol);
nurbs = cell(1,noPatches);


for patch = 1:noPatches
    nurbs_patch = cell(1,n);
    for i = 1:n
        nurbs_patch(i) = nurbsCol{i}(patch);
    end
    nurbs_patch = homogenizeNURBSparametrization(nurbs_patch); % Ensure that the basis is the same

    coeffs = zeros([size(nurbs_patch{1}.coeffs), n]);
    for i = 1:n
        coeffs_i = nurbs_patch{i}.coeffs;
        subs = {[repmat({':'},1,ndims(coeffs)-1),{i}]};
        coeffs = subasgnArr(coeffs,coeffs_i,subs);
    end
    Eta = 0:(n-p);
    Eta = [zeros(1,p),Eta/(n-p),ones(1,p)];
    nurbs(patch) = createNURBSobject(coeffs,[nurbs_patch{1}.knots(:)', {Eta}]);
end
% Change the parametric direction of lofting to be aligned with the parametric direction dir
switch dir
    case 1
        nurbs = permuteNURBS(nurbs,[3,1,2]);
    case 2
        nurbs = permuteNURBS(nurbs,[2,3,1]);
end