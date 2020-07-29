function nurbs = scaleNURBSweights(nurbs,scale)

if nargin < 2
    scale = -Inf;
    for patch = 1:numel(nurbs)
        w = max(nurbs{patch}.coeffs(4,:));
        if w > scale
            scale = w;
        end
    end
end
for patch = 1:numel(nurbs)
    nurbs{patch}.coeffs(4,:) = nurbs{patch}.coeffs(4,:)/scale;
end
