function X = computeBoundingBox(nurbs)

X = zeros(3,2);
X(:,1) = Inf;
X(:,2) = -Inf;
for i = 1:numel(nurbs)
    for j = 1:nurbs{i}.d
        if X(j,1) > min(nurbs{i}.coeffs(j,:))
            X(j,1) = min(nurbs{i}.coeffs(j,:));
        end
        if X(j,2) < max(nurbs{i}.coeffs(j,:))
            X(j,2) = max(nurbs{i}.coeffs(j,:));
        end
    end
end
    