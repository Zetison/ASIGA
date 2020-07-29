function fv = fixBeTSSiTriangularization(fv,alpha_threshold)

c = 4;
s = 1.2; 
Eps = 1e-6;
P = fv.vertices;
% [~, gluedNodes] = uniquetol(P,Eps,'ByRows',true, 'DataScale',max(norm2(P)), 'OutputAllIndices', true);
% [P, IA,IC] = uniquetol(P,Eps,'ByRows',true, 'DataScale',max(norm2(P)));
% tri = 
tri = fv.faces;
% noElems = size(tri,1);
% dofs = size(P,1);

l1 = norm2(P(tri(:,2),:)-P(tri(:,3),:));
l2 = norm2(P(tri(:,1),:)-P(tri(:,3),:));
l3 = norm2(P(tri(:,1),:)-P(tri(:,2),:));
% h_max = max(2*l1.*l2.*l3./sqrt((l1+l2+l3).*(l1+l2-l3).*(l1+l3-l2).*(l2+l3-l1)));
% l_min = min([l1, l2, l3],[],2);
% l_max = max([l1, l2, l3],[],2);
% aspect = max(l_max./l_min);

alpha1 = 180*acos((l2.^2+l3.^2-l1.^2)./(2*l2.*l3))/pi;
alpha2 = 180*acos((l1.^2+l3.^2-l2.^2)./(2*l1.*l3))/pi;
alpha3 = 180*acos((l1.^2+l2.^2-l3.^2)./(2*l1.*l2))/pi;
alpha = [alpha1, alpha2, alpha3];
badtri = find(max(alpha,[],2) > alpha_threshold);
z = [P(tri(badtri,1),3),P(tri(badtri,2),3),P(tri(badtri,3),3)];
i2 = find(any(abs(z-c) < Eps,2));
y = [P(tri(badtri(i2),1),2),P(tri(badtri(i2),2),2),P(tri(badtri(i2),3),2)];
badtri2 = badtri(i2(any(abs(abs(y)-s) < Eps,2)));
[~, ind] = min(alpha(badtri2,:), [], 2);
for i = 1:numel(badtri2)
    indices = setdiff([1,2,3],ind(i));
    idx1 = tri(badtri2(i),indices(1));
    idx2 = tri(badtri2(i),indices(2));
    if abs(P(idx1,3) - c) < Eps && abs(abs(P(idx1,2)) - s) < Eps
        P(idx1,1) = P(idx2,1);
    elseif abs(P(idx2,3) - c) < Eps && abs(abs(P(idx2,2)) - s) < Eps
        P(idx2,1) = P(idx1,1);
    end
end
fv.vertices = P;
% P(tri(indices

% f = 1000;
% omega = 2*pi*f;
% k = omega/1500;
% lambda = 2*pi/k;
end
