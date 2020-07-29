function [noElems, dofs, h_max1, h_max2, alpha_max, alpha_min, aspect, waterTight, area_min, area_max, skewness, i, Pi] = computeMeshData(filename)

if isstring(filename) || isa(filename,'char')
    fv = stlread(filename); 
else
    fv = filename;
end
tri = fv.faces;
P = fv.vertices;
noElems = size(tri,1);
dofs = size(P,1);
P23 = P(tri(:,2),:)-P(tri(:,3),:);
P13 = P(tri(:,1),:)-P(tri(:,3),:);
l1 = norm2(P23);
l2 = norm2(P13);
l3 = norm2(P(tri(:,1),:)-P(tri(:,2),:));
l_min = min([l1, l2, l3],[],2);
l_max = max([l1, l2, l3],[],2);
aspect = l_max./l_min;
h_max1 = max(2*l1.*l2.*l3./sqrt((l1+l2+l3).*(l1+l2-l3).*(l1+l3-l2).*(l2+l3-l1)));
h_max2 = max(l_max);
% h_max = max((l1 + l2 + l3)/3);

areas = norm2(cross(P23,P13,2))/2;
area_min = min(areas);
area_max = max(areas);
if 0
    idx1 = tri(:,1);
    idx2 = tri(:,2);
    waterTight = true(noElems,1);
    parfor i = 1:noElems
        if ~(numel(find(and(any(tri==idx1(i),2),any(tri==idx2(i),2)))) == 2)
            waterTight(i) = false;
        end
    end
    waterTight = ~any(~waterTight);
    if ~waterTight
        warning('Mesh is not water tight!')
    end
else
    waterTight = true;
end
if area_min < eps
    warning('Mesh contains triangles of measure zero!')
end
alpha1 = 180*acos((l2.^2+l3.^2-l1.^2)./(2*l2.*l3))/pi;
alpha2 = 180*acos((l1.^2+l3.^2-l2.^2)./(2*l1.*l3))/pi;
alpha3 = 180*acos((l1.^2+l2.^2-l3.^2)./(2*l1.*l2))/pi;
alpha = [alpha1, alpha2, alpha3];
[alpha_max, i] = max(max(alpha,[],2));
alpha_min = min(min(alpha,[],2));
Pi = P(tri(i,:),:);

alpha_e = 60;
alpha = alpha(:);
skewness = min(1-max([(alpha-alpha_e)/(180-alpha_e), (alpha_e-alpha)/alpha_e],[],2));
% f = 1000;
% omega = 2*pi*f;
% k = omega/1500;
% lambda = 2*pi/k;
end
