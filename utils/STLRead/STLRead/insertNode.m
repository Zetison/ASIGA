function fv = insertNode(fv,alpha_threshold)

P = fv.vertices;
tri = fv.faces;
noElems = size(tri,1);

l1 = norm2(P(tri(:,2),:)-P(tri(:,3),:));
l2 = norm2(P(tri(:,1),:)-P(tri(:,3),:));
l3 = norm2(P(tri(:,1),:)-P(tri(:,2),:));

alpha1 = 180*acos((l2.^2+l3.^2-l1.^2)./(2*l2.*l3))/pi;
alpha2 = 180*acos((l1.^2+l3.^2-l2.^2)./(2*l1.*l3))/pi;
alpha3 = 180*acos((l1.^2+l2.^2-l3.^2)./(2*l1.*l2))/pi;
alpha = [alpha1, alpha2, alpha3];
badtri = find(max(alpha,[],2) > alpha_threshold);
[~, ind] = max(alpha(badtri,:), [], 2);
no_badtri = numel(badtri);
tri2 = [tri; zeros(4*no_badtri,3)];
P2 = [P; zeros(no_badtri,3)];
badtri2 = zeros(no_badtri,1);
for i = 1:numel(badtri)
    indices = setdiff([1,2,3],ind(i));
    i_badtri = badtri(i);
    idx1 = tri(i_badtri,indices(1));
    idx2 = tri(i_badtri,indices(2));
    badtri2(i) = setdiff(find(and(any(tri==idx1,2),any(tri==idx2,2))),badtri(i));
end
counter = noElems+1;
counter2 = 1;
badtri3 = zeros(no_badtri,2);
for i = 1:numel(badtri)
    i_badtri2 = badtri2(i);
    if ~ismember(i_badtri2,badtri2(i+1:end)) % no other bad triangles claim to divide the neighbouring triangle
        indices = setdiff([1,2,3],ind(i));
        i_badtri = badtri(i);
        idx1 = tri2(i_badtri,indices(1));
        idx2 = tri2(i_badtri,indices(2));
        P_new = (P2(idx1,:)+P2(idx2,:))/2;
        P2(idx1+1:end,:) = P2(idx1:end-1,:);
        P2(idx1,:) = P_new;
        indicesTri2 = and(tri2 >= idx1, tri2 ~= 0);
        tri2(indicesTri2) = tri2(indicesTri2) + 1;
        
        temp = tri2(i_badtri,:);
        temp(indices(1)) = idx1;
        tri2(counter,:) = temp;
        temp = tri2(i_badtri,:);
        temp(indices(2)) = idx1;
        tri2(counter+1,:) = temp;
        
        temp = tri2(i_badtri2,:);
        temp(temp == tri2(i_badtri,indices(1))) = idx1;
        tri2(counter+2,:) = temp;
        temp = tri2(i_badtri2,:);
        temp(temp == tri2(i_badtri,indices(2))) = idx1;
        tri2(counter+3,:) = temp;
        
        badtri3(counter2,:) = [i_badtri,i_badtri2];
        counter2 = counter2 + 1;
        counter = counter + 4;
    end
end
fprintf('Inserted %d new nodes. There were %d bad triangles\n', counter2-1,no_badtri);
badtri3(counter2:end,:) = [];
tri2(counter:end,:) = [];
tri2(unique(badtri3(:)),:) = [];
fv.vertices = P2;
fv.faces = sortrows(tri2,1);

end
