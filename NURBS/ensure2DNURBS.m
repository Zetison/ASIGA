function nurbs = ensure2DNURBS(nurbs)
if ~iscell(nurbs)
    nurbs = {nurbs};
end
for i = 1:numel(nurbs)
    coeffs = nurbs{i}.coeffs;
    d = nurbs{i}.d;
    if d == 1
        sizes = size(coeffs);
        nurbs{i}.coeffs = cat(1,slc(coeffs,1:d),zeros([2-d,sizes(2:end)]),slc(coeffs,d+1));
        nurbs{i}.d = 2;
    elseif d > 2
        nurbs{i}.coeffs = slc(nurbs{i}.coeffs, [1:2,size(nurbs{i}.coeffs,1)], 1);
        nurbs{i}.d = 2;
    end
end
