function patchTop = getPatchTopology(patches)
% This routine assumes all patches have same orientation in the sense that
% the normal vector points in the same direction
noPatches = numel(patches);
patchData = cell(noPatches,1);
Eps = 1e7*eps;
for i = 1:noPatches
    patchData{i}.edgeNodes = cell(4,1);
    patchData{i}.edgeNodes{1} = permute(patches{i}.coeffs(:,end,:),[1,3,2]);
    patchData{i}.edgeNodes{2} = patches{i}.coeffs(:,:,end);
    patchData{i}.edgeNodes{3} = permute(patches{i}.coeffs(:,1,:),[1,3,2]);
    patchData{i}.edgeNodes{4} = patches{i}.coeffs(:,:,1);
end
patchTop = cell(noPatches,1);
for i = 1:noPatches
    patchTop{i} = NaN(4,2);
    for dir_i = 0:3
        pts_i = patchData{i}.edgeNodes{dir_i+1};
        if size(uniquetol(pts_i(1:3,:).',Eps,'ByRows',true, 'DataScale',max(norm2(pts_i(1:3,:).'))),1) == 1 % singular point
            continue
        end
        foundSide = false;
        for j = 1:noPatches
            for dir_j = 0:3
                pts_j = patchData{j}.edgeNodes{dir_j+1};
                if size(pts_i,2) ~= size(pts_j,2)
                    continue
                end
                flip_j_side = false;
                twist = mod(dir_i-dir_j+2,4);
                if mod(dir_i,2)
                    if twist == 1 || twist == 2
                        flip_j_side = true;
                    end
                else
                    if twist == 2 || twist == 3
                        flip_j_side = true;
                    end
                end
                if flip_j_side
                    diff = abs(pts_i-fliplr(pts_j));
                else
                    diff = abs(pts_i-pts_j);
                end
                if all(diff(:) < Eps) && ~(i == j && dir_i == dir_j)
                    patchTop{i}(dir_i+1,:) = [j, twist];
                    foundSide = true;
                    break
                end
            end
            if foundSide
                break
            end
        end
    end
end