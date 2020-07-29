function nurbs = loftNURBS(nurbsCol)
% This function only give a linear lofting between curves
for i = 1:numel(nurbsCol)
    nurbsCol{i} = ensure3DNURBS(nurbsCol{i});
end
noPatches = numel(nurbsCol{1});
n = numel(nurbsCol);
nurbs = cell(1,noPatches);
for patch = 1:noPatches
    coeffs = zeros([size(nurbsCol{1}{patch}.coeffs), n]);
    for i = 1:n
        coeffs_i = nurbsCol{i}{patch}.coeffs;
        subs = {[repmat({':'},1,ndims(coeffs)-1),{i}]};
        coeffs = subasgnArr(coeffs,coeffs_i,subs);
    end
    Eta = 0:(n-1);
    Eta = [0,Eta/(n-1),1];
    nurbs(patch) = createNURBSobject(coeffs,[nurbsCol{1}{patch}.knots(:)', {Eta}]);
end