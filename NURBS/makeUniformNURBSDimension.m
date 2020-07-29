function nurbs = makeUniformNURBSDimension(nurbs, d)

d_max = -Inf;
for i = 1:numel(nurbs)
    if nurbs{i}.d > d_max
        d_max = nurbs{i}.d;
    end
end
if nargin > 1
    d_max = max(d_max,d);
end
for i = 1:numel(nurbs)
    d = size(nurbs{i}.coeffs,1)-1;
    nurbs{i}.coeffs(d_max+1,:) = nurbs{i}.coeffs(d+1,:);
    nurbs{i}.coeffs(d+1:d_max,:) = 0;
end