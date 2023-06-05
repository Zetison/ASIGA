function nurbs_vol = surfaceToVolume(nurbs_surf,connection,t)
% This routine assumes nurbs_surf to be a closed NURBS-surface with all
% patches having normal vectors pointing outwards.
d_p = nurbs_surf{1}.d_p;
d = 3; % dimension

noBdryPatches = numel(nurbs_surf);

topologyMap = createTopologyMap(connection,noBdryPatches,d_p);
[cornerData, angles] = computeCornerData(nurbs_surf,topologyMap);

noNakedNurbsPatches = noBdryPatches;
nurbs_bdry = nurbs_surf;
idxTo_idx1_orient = [0,3,5,6];
idxTo_idx3_orient = [6,5,0,3];
idxTo_idx3_orient_flipped = [4,7,1,2];

nurbs_vol = cell(1,4*numel(nurbs_surf));
counter = 1;
midx2corners = [1,4; 2,3; 1,2; 4,3];
% sharpAngle = 120*pi/180; % Threshold for a "sharp" angle
sharpAngle = 125*pi/180; % Threshold for a "sharp" angle
avg_v_n_threshholdAngle = 11.25*pi/180; % Threshold angle deviation between the normal vectors v_n avg_v_n
original_angles = angles;

while any(angles(:) < pi) && noNakedNurbsPatches
    [patch,maxNoSharpAngles,minSum] = findNextPatch(original_angles,sharpAngle); % Check for sharp angles in at the initial nurbs_surf first
    midx = find(original_angles(patch,:) < sharpAngle);
    if isinf(minSum)
        [patch,maxNoSharpAngles] = findNextPatch(angles,sharpAngle);
        midx = find(angles(patch,:) < sharpAngle);
    end
    faces = cell(1,6);

%     for i = 1:numel(topologyMap(patch))
    switch maxNoSharpAngles
        case 0
            g = cell(1,d_p);
            for i = 1:numel(knots)
                g{i} = aveknt(knots{i}, degree(i)+1);
            end
            
            X = evaluateNURBSvec(nurbs_bdry,[g]);
            nurbs_covered = patch;
        case 1 % Loft path along slave
            slave = topologyMap{patch}(midx).slave;
            sidx = topologyMap{patch}(midx).sidx;
            faces(1) = nurbs_bdry(patch);
            faces(3) = nurbs_bdry(slave);

            faces{1} = orientNURBS(faces{1}, idxTo_idx1_orient(midx)); % Make "midx = 1"
            faces{3}  = orientNURBS(faces{3},  idxTo_idx3_orient(sidx));   % Make "sidx = 3"
            knots = [faces{3}.knots(end), faces{1}.knots];
            degree = [faces{3}.degree(end), faces{1}.degree];
            number = [faces{3}.number(end), faces{1}.number];
            master_coeffs = faces{1}.coeffs;
            master_coeffs = reshape(master_coeffs,[d+1,1,faces{1}.number]);
            slave1_coeffs = permute(faces{3}.coeffs,[1,3,2]);
            slave1_coeffs = reshape(slave1_coeffs,[d+1,faces{3}.number(2),1,faces{3}.number(1)]);
            
            midx_opposite = midx-(-1)^midx;
            if all(angles(patch,midx) >= pi/2)
%             if angles(patch,midx) + angles(patch,midx_opposite)/2 > pi % Make patching normal based (less fear of intersection) % TODO: REIMPLEMENT THIS WITH TRY AND CATCH (for when Jacobian is negative - intersection)
                oppositeCornerIndices = midx2corners(midx_opposite,:);
                if midx_opposite == 1 || midx_opposite == 4
                    oppositeCornerIndices = fliplr(oppositeCornerIndices);
                end
                cornerNormal = cornerData(patch).avg_v_n(oppositeCornerIndices,:).';
                v_n = cornerData(patch).v_n(oppositeCornerIndices,:).';
                indices = acos(abs(dot(cornerNormal,v_n,1))) > avg_v_n_threshholdAngle;
                cornerNormal(:, indices) = v_n(:,indices);
                len1 = norm(slave1_coeffs(1:3,end,1,1) - slave1_coeffs(1:3,1,1,1));
                len2 = norm(slave1_coeffs(1:3,end,1,end) - slave1_coeffs(1:3,1,1,end));
                g = reshape(aveknt(knots{1}, degree(1)+1),1,[]);
                coeffs1 = zeros(d+1,number(1));
                coeffs1(1:3,:) = master_coeffs(1:3,1,end,1) + len1*cornerNormal(:,1).*g;
                coeffs1(4,:) = [master_coeffs(4,1,end,1), slave1_coeffs(4,2:end,1,1)]; % Use the weights based on the weights on opposite side
                coeffs2 = zeros(d+1,number(1));
                coeffs2(1:3,:) = master_coeffs(1:3,1,end,end) + len2*cornerNormal(:,2).*g;
                coeffs2(4,:) = [master_coeffs(4,1,end,end), slave1_coeffs(4,2:end,1,end)]; % Use the weights based on the weights on opposite side
                g = reshape(aveknt(knots{3}, degree(3)+1),1,[]);
                coeffs4 = coeffs1(:,end).*(1-g) + coeffs2(:,end).*g;
    
                edges = cell(1,4);
                edges(1) = createNURBSobject(coeffs1,knots{1});
                edges(2) = createNURBSobject(coeffs2,knots{1});
                edges(3) = subNURBS(faces(1),'at',[0,1;0,0]);
                edges(4) = createNURBSobject(coeffs4,knots{3});
                faces(4) = GordonHall(edges);
                edges_s = subNURBS(faces(3),'at',[0,0;0,1]);
                faces(2) = loftNURBS({edges{4},edges_s{1}});
                faces(2) = elevateNURBSdegree(faces(2),[0,1]);
                faces(2) = permuteNURBS(faces(2),[2,1]);
                edges(1) = subNURBS(faces(1),'at',[0,0;1,0]);
                edges(2) = subNURBS(faces(2),'at',[0,0;1,0]);
                edges(3) = subNURBS(faces(3),'at',[1,0;0,0]);
                edges(4) = subNURBS(faces(4),'at',[1,0;0,0]);
                faces(5) = GordonHall(edges);
                edges(1) = subNURBS(faces(1),'at',[0,0;0,1]);
                edges(2) = subNURBS(faces(2),'at',[0,0;0,1]);
                edges(3) = subNURBS(faces(3),'at',[0,1;0,0]);
                edges(4) = subNURBS(faces(4),'at',[0,1;0,0]);
                faces(6) = GordonHall(edges);
                
                nurbs_vol(counter) = GordonHall(faces);
                
%                 if checkOrientation(nurbs_vol(counter), 10)
            else
                coeffs = master_coeffs + slave1_coeffs - slave1_coeffs(:,1,1,:);
                nurbs_vol(counter) = createNURBSobject(coeffs,knots);
            end
            nurbs_covered = [patch,slave];
            nurbs_newBdry = subNURBS(nurbs_vol(counter),'at',[0,1;0,1;1,1],'outwardPointingNormals',true);
        case 2 % Loft path along slave
            if all(ceil(midx/2) == [1,2]) % Fill a patch in a corner
            
                if all(midx == [2,3])
                    midx = [3,2];
                elseif all(midx == [1,4])
                    midx = [4,1];
                end
                slave1 = topologyMap{patch}(midx(1)).slave;
                sidx1 = topologyMap{patch}(midx(1)).sidx;
                slave2 = topologyMap{patch}(midx(2)).slave;
                sidx2 = topologyMap{patch}(midx(2)).sidx;

                faces(1) = nurbs_bdry(patch);
                faces(3) = nurbs_bdry(slave1);
                faces(5) = nurbs_bdry(slave2);

                faces{1} = orientNURBS(faces{1}, idxTo_idx1_orient(midx(1))); % Make "midx = 1"
                faces{3} = orientNURBS(faces{3}, idxTo_idx3_orient(sidx1));   % Make "sidx = 3"
                faces{5} = orientNURBS(faces{5}, idxTo_idx1_orient(sidx2));   % Make "sidx = 1"


                knots = [faces{3}.knots(end), faces{1}.knots];
                degree = [faces{3}.degree(end), faces{1}.degree];
                number = [faces{3}.number(end), faces{1}.number];
                master_coeffs = faces{1}.coeffs;
                master_coeffs = reshape(master_coeffs,[d+1,1,faces{1}.number]);
                slave1_coeffs = permute(faces{3}.coeffs,[1,3,2]);
                slave1_coeffs = reshape(slave1_coeffs,[d+1,faces{3}.number(2),1,faces{3}.number(1)]);
                slave2_coeffs = faces{5}.coeffs;
                midx_opposite = midx-(-1).^midx;
%                 if all(angles(patch,midx) + angles(patch,midx_opposite)/2 > pi) % Make patching normal based (less fear of intersection) % TODO: REIMPLEMENT THIS WITH TRY AND CATCH (for when Jacobian is negative - intersection)
                if all(angles(patch,midx) >= pi/2) % Make patching normal based (less fear of intersection) % TODO: REIMPLEMENT THIS WITH TRY AND CATCH (for when Jacobian is negative - intersection)
                    oppositeCornerIndices = midx2corners(midx_opposite(1),:);
                    if midx_opposite(1) == 1 || midx_opposite(1) == 4
                        oppositeCornerIndices = fliplr(oppositeCornerIndices);
                    end
                    cornerNormal = cornerData(patch).avg_v_n(oppositeCornerIndices,:).';
                    v_n = cornerData(patch).v_n(oppositeCornerIndices,:).';
                    indices = acos(abs(dot(cornerNormal,v_n,1))) > avg_v_n_threshholdAngle;
                    cornerNormal(:, indices) = v_n(:,indices);
                    len2 = norm(slave1_coeffs(1:3,end,1,end) - slave1_coeffs(1:3,1,1,end));
                    g = reshape(aveknt(knots{1}, degree(1)+1),1,[]);
                    coeffs2 = zeros(d+1,number(1));
                    coeffs2(1:3,:) = master_coeffs(1:3,1,end,end) + len2*cornerNormal(:,2).*g;
                    coeffs2(4,:) = [master_coeffs(4,1,end,end), slave1_coeffs(4,2:end,1,end)]; % Use the weights based on the weights on opposite side
    
                    edges = cell(1,4);
                    edges(1) = subNURBS(faces(5),'at',[0,0;0,1]);
                    edges(2) = createNURBSobject(coeffs2,knots{1});
                    edges(3) = subNURBS(faces(1),'at',[0,1;0,0]);
                    g = reshape(aveknt(knots{3}, degree(3)+1),1,[]);
                    coeffs4 = edges{1}.coeffs(:,end).*(1-g) + coeffs2(:,end).*g;
                    edges(4) = createNURBSobject(coeffs4,knots{3});
                    faces(4) = GordonHall(edges);
    
                    edges(1) = subNURBS(faces(3),'at',[0,0;0,1]);
                    edges(2) = subNURBS(faces(4),'at',[0,0;0,1]);
                    edges(3) = subNURBS(faces(5),'at',[0,1;0,0]);
                    g = reshape(aveknt(knots{2}, degree(2)+1),1,[]);
                    coeffs4 = edges{1}.coeffs(:,end).*(1-g) + edges{2}.coeffs(:,end).*g;
                    edges(4) = createNURBSobject(coeffs4,knots{2});
                    faces(2) = GordonHall(edges);
    
                    edges(1) = subNURBS(faces(1),'at',[0,0;0,1]);
                    edges(2) = subNURBS(faces(2),'at',[0,0;0,1]);
                    edges(3) = subNURBS(faces(3),'at',[0,1;0,0]);
                    edges(4) = subNURBS(faces(4),'at',[0,1;0,0]);
                    faces(6) = GordonHall(edges);
                    
                    nurbs_vol(counter) = GordonHall(faces);

%                 if checkOrientation(nurbs_vol(counter), 10)
                else
                    coeffs = master_coeffs + slave1_coeffs + slave2_coeffs - master_coeffs(:,1,1,:) - slave2_coeffs(:,1,:,1) - slave1_coeffs(:,:,1,1) + master_coeffs(:,1,1,1);
        
                    nurbs_vol(counter) = createNURBSobject(coeffs,knots);
                end
                nurbs_covered = [patch,slave1,slave2];
                nurbs_newBdry = subNURBS(nurbs_vol(counter),'at',[0,1;0,1;0,1],'outwardPointingNormals',true);
            else
                slave1 = topologyMap{patch}(midx(1)).slave;
                sidx1 = topologyMap{patch}(midx(1)).sidx;
                slave2 = topologyMap{patch}(midx(2)).slave;
                sidx2 = topologyMap{patch}(midx(2)).sidx;

                faces(1) = nurbs_bdry(patch);
                faces(3) = nurbs_bdry(slave1);
                faces(4) = nurbs_bdry(slave2);

                faces{1} = orientNURBS(faces{1}, idxTo_idx1_orient(midx(1))); % Make "midx = 1"
                faces{3} = orientNURBS(faces{3}, idxTo_idx3_orient(sidx1));   % Make "sidx = 3"
                faces{4} = orientNURBS(faces{4}, idxTo_idx3_orient_flipped(sidx2));   % Make "sidx = 1"


                knots = [faces{3}.knots(end), faces{1}.knots];
                master_coeffs = faces{1}.coeffs;
                master_coeffs = reshape(master_coeffs,[d+1,1,faces{1}.number]);
                slave1_coeffs = permute(faces{3}.coeffs,[1,3,2]);
                slave1_coeffs = reshape(slave1_coeffs,[d+1,faces{3}.number(2),1,faces{3}.number(1)]);
                slave2_coeffs = permute(faces{4}.coeffs,[1,3,2]);
                slave2_coeffs = reshape(slave2_coeffs,[d+1,faces{4}.number(2),1,faces{4}.number(1)]);
                g = reshape(aveknt(faces{1}.knots{1},faces{1}.degree(1)+1),1,1,[]);
                coeffs = (master_coeffs + slave1_coeffs - slave1_coeffs(:,1,1,:)).*(1-g) + (master_coeffs + slave2_coeffs - slave2_coeffs(:,1,1,:)).*g;
    
                nurbs_vol(counter) = createNURBSobject(coeffs,knots);

                nurbs_covered = [patch,slave1,slave2];
                nurbs_newBdry = subNURBS(nurbs_vol(counter),'at',[0,1;0,0;1,1],'outwardPointingNormals',true);
            end
        case 3
            if all(midx == [1,2,3])
                midx = [1,3,2];
            elseif all(midx == [1,2,4])
                midx = [2,4,1];
            elseif all(midx == [1,3,4])
                midx = [4,1,3];
            elseif all(midx == [2,3,4])
                midx = [3,2,4];
            end
            slave1 = topologyMap{patch}(midx(1)).slave;
            sidx1 = topologyMap{patch}(midx(1)).sidx;
            slave2 = topologyMap{patch}(midx(2)).slave;
            sidx2 = topologyMap{patch}(midx(2)).sidx;
            slave3 = topologyMap{patch}(midx(3)).slave;
            sidx3 = topologyMap{patch}(midx(3)).sidx;

            faces(1) = nurbs_bdry(patch);
            faces(3) = nurbs_bdry(slave1);
            faces(5) = nurbs_bdry(slave2);
            faces(4) = nurbs_bdry(slave3);

            faces{1} = orientNURBS(faces{1}, idxTo_idx1_orient(midx(1))); % Make "midx = 1"
            faces{3} = orientNURBS(faces{3}, idxTo_idx3_orient(sidx1));   % Make "sidx = 3"
            faces{5} = orientNURBS(faces{5}, idxTo_idx1_orient(sidx2));   % Make "sidx = 1"
            faces{4} = orientNURBS(faces{4}, idxTo_idx3_orient_flipped(sidx3));   % Make "sidx = 1"


            knots = [faces{3}.knots(end), faces{1}.knots];
            master_coeffs = faces{1}.coeffs;
            master_coeffs = reshape(master_coeffs,[d+1,1,faces{1}.number]);
            slave1_coeffs = permute(faces{3}.coeffs,[1,3,2]);
            slave1_coeffs = reshape(slave1_coeffs,[d+1,faces{3}.number(2),1,faces{3}.number(1)]);
            slave2_coeffs = faces{5}.coeffs;
            slave3_coeffs = permute(faces{4}.coeffs,[1,3,2]);
            slave3_coeffs = reshape(slave3_coeffs,[d+1,faces{4}.number(2),1,faces{4}.number(1)]);
            g = reshape(aveknt(faces{1}.knots{1},faces{1}.degree(1)+1),1,1,[]);
            coeffs = (master_coeffs + slave1_coeffs + slave2_coeffs - master_coeffs(:,1,1,:) - slave2_coeffs(:,1,:,1) - slave1_coeffs(:,:,1,1) + master_coeffs(:,1,1,1)).*(1-g) ...
                   + (master_coeffs + slave3_coeffs + slave2_coeffs - master_coeffs(:,1,end,:) - slave2_coeffs(:,1,:,1) - slave3_coeffs(:,:,end,1) + master_coeffs(:,1,end,1)).*g;

            nurbs_vol(counter) = createNURBSobject(coeffs,knots);

            nurbs_covered = [patch,slave1,slave2,slave3];
            nurbs_newBdry = subNURBS(nurbs_vol(counter),'at',[0,1;0,0;0,1],'outwardPointingNormals',true);
        case 4
            nurbs_vol(counter) = GordonHall(faces);
        otherwise
            error('not implemented')
    
    end

    % Remove faces with zero measure
    zeroMeasure = NURBShasZeroMeasure(nurbs_newBdry);
    nurbs_newBdry(zeroMeasure) = [];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Plotting
    plotNURBSvec(nurbs_vol(counter),'plotControlPolygon',0,'plotNormalVectors',0);
    drawnow
    pause(0.001)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    noNewBdryPatches = numel(nurbs_newBdry);
    newBdryPatches = (numel(nurbs_bdry)+1):(numel(nurbs_bdry)+noNewBdryPatches);

    % Update nurbs_bdry
    nurbs_bdry = [nurbs_bdry, nurbs_newBdry];
    surroundingPatches = 4*numel(nurbs_covered);
    counter2 = 1;
    for i = 1:numel(nurbs_covered)
        patch_i = nurbs_covered(i);
        for midx_i = 1:numel(topologyMap{patch_i})
            if ~isempty(topologyMap{patch_i}(midx_i).slave)
                surroundingPatches(counter2) = topologyMap{patch_i}(midx_i).slave;
                counter2 = counter2 + 1;
            end
        end
    end
    surroundingPatches(counter2:end) = [];
    surroundingPatches = setdiff(surroundingPatches,nurbs_covered); % Remove the patches that is replaced
    local2global = [surroundingPatches,newBdryPatches];
    geometry = getTopology(nurbs_bdry(local2global));
    connectionSub = geometry.topology.connection;

    for patch_i = newBdryPatches
        topologyMap{patch_i} = struct('slave', cell(1, 4), ...
                                      'sidx', cell(1, 4), ...
                                      'orient', cell(1, 4));
    end

    % Update topologyMap
    for i = 1:numel(connectionSub)
        patch_i = local2global(connectionSub{i}.Attributes.master);
        slave_i = local2global(connectionSub{i}.Attributes.slave);
        midx_i = connectionSub{i}.Attributes.midx;
        sidx_i = connectionSub{i}.Attributes.sidx;
        orient = connectionSub{i}.Attributes.orient;

        topologyMap{patch_i}(midx_i).slave = slave_i;
        topologyMap{patch_i}(midx_i).sidx = sidx_i;
        topologyMap{patch_i}(midx_i).orient = orient;
    end

    % Update angles
    [cornerData,angles] = computeCornerData(nurbs_bdry,topologyMap,newBdryPatches,cornerData,angles);
    angles(nurbs_covered,:) = NaN;

    noNakedNurbsPatches = noNakedNurbsPatches - 1 - maxNoSharpAngles;
    counter = counter + 1;
    if patch <= noBdryPatches
        original_angles(patch,:) = NaN;
        for i_midx = 1:numel(midx)
            slave = topologyMap{patch}(midx(i_midx)).slave;
            original_angles(slave,:) = NaN;
        end
    end
end
nurbs_vol(counter:end) = [];

function [patch,maxNoSharpAngles,minSum] = findNextPatch(angles,sharpAngle)

sharpAngles = angles < sharpAngle;
noSharpAngles = sum(sharpAngles,2);
maxNoSharpAngles = max(noSharpAngles);
candidates = noSharpAngles == maxNoSharpAngles;
avgAngles = Inf(size(sharpAngles,1),1);
for i = 1:size(sharpAngles,1)
    if candidates(i)
        angles_i = angles(i,:);
        angles_i(isnan(angles_i)) = [];
%             avgAngles(i) = sum(angles_i)/numel(angles_i);
        avgAngles(i) = min(angles_i);
    end
end

[minSum,patch] = min(avgAngles);

function [cornerData,angles] = computeCornerData(nurbs,topologyMap,indices,cornerData,angles)
noBdryPatches = numel(nurbs);
d_p = nurbs{1}.d_p;
d = 3; % dimension
if nargin < 3
    indices = 1:noBdryPatches;
end
if nargin < 4
    cornerData = struct('X', cell(1, noBdryPatches), ...
                        'v_t1', cell(1, noBdryPatches), ...
                        'v_t2', cell(1, noBdryPatches), ...
                        'v_n', cell(1, noBdryPatches), ...
                        'avg_v_n', cell(1, noBdryPatches));
    angles = NaN(noBdryPatches,2*d_p);
else
    noNewFields = noBdryPatches-numel(cornerData);
    cornerData = [cornerData, struct('X', cell(1, noNewFields), ...
                                    'v_t1', cell(1, noNewFields), ...
                                    'v_t2', cell(1, noNewFields), ...
                                    'v_n', cell(1, noNewFields), ...
                                    'avg_v_n', cell(1, noNewFields))];
    angles = [angles; NaN(noNewFields,2*d_p)];
end
Eps = 1e-10;


xi = [0,0;
      1,0;
      1,1;
      0,1]; % Corner indices are oriented counter-clockwise starting at xi,eta=0
xi_l = 0.75; % a somewhat arbitrary value different from 0,0.5,1 (to avoid singular problems)

%% Compute normal vector at the vertices
for patch = indices
    [X,dXdxi,dXdeta] = evaluateNURBSvec(nurbs{patch},xi,1);
    cornerData(patch).X = X;
    v_t1 = dXdxi./norm2(dXdxi);
    v_t2 = dXdeta./norm2(dXdeta);
    cornerData(patch).v_t1 = v_t1;
    cornerData(patch).v_t2 = v_t2;
    if norm(X(1,:)-X(2,:)) < Eps
        [~,~,dXdeta_m] = evaluateNURBSvec(nurbs{patch},[xi_l,0],1);
        v_t2_m = dXdeta_m./norm(dXdeta_m);
        v_t1(1,:) = v_t2_m;
        v_t1(2,:) = -v_t2_m;
    end
    if norm(X(2,:)-X(3,:)) < Eps
        [~,dXdxi_m,~] = evaluateNURBSvec(nurbs{patch},[1,xi_l],1);
        v_t1_m = dXdxi_m./norm(dXdxi_m);
        v_t2(2,:) = -v_t1_m;
        v_t2(3,:) = v_t1_m;
    end
    if norm(X(3,:)-X(4,:)) < Eps
        [~,~,dXdeta_m] = evaluateNURBSvec(nurbs{patch},[xi_l,1],1);
        v_t2_m = dXdeta_m./norm(dXdeta_m);
        v_t1(3,:) = v_t2_m;
        v_t1(4,:) = -v_t2_m;
    end
    if norm(X(4,:)-X(1,:)) < Eps
        [~,dXdxi_m,~] = evaluateNURBSvec(nurbs{patch},[0,xi_l],1);
        v_t1_m = dXdxi_m./norm(dXdxi_m);
        v_t2(4,:) = -v_t1_m;
        v_t2(1,:) = v_t1_m;
    end
    cornerData(patch).v_n = NaN(4,d);
    v_n = cross(v_t1, v_t2, 2);
    vectorLengths = norm2(v_n);
    indices_Lipschitz = vectorLengths > Eps; % skip the cases of "non-Lipschitz vertices"
    cornerData(patch).v_n(indices_Lipschitz,:) = v_n(indices_Lipschitz,:)./vectorLengths(indices_Lipschitz);
end
% for patch = 1:noBdryPatches
%     quiver3(cornerData(patch).X(:,1), ...
%             cornerData(patch).X(:,2), ...
%             cornerData(patch).X(:,3), ...
%             cornerData(patch).v_n(:,1), ...
%             cornerData(patch).v_n(:,2), ...
%             cornerData(patch).v_n(:,3),'LineWidth',1,'AutoScale','off','DisplayName',['Patch ' num2str(patch)])
% end

%% Compute average normal vector in all vertices
noCorners = 2^d_p;
corner2midx = [1,3,2,4]; % This map gives the midx for the side of the patch sharing a given corner. 
sidx2corner = [4,2,1,3]; % This map gives the i_corner for the corner of the patch from a given side. 
for patch = indices
    cornerData(patch).avg_v_n = zeros(size(v_t1));
    for i_corner = 1:noCorners
        v_n = zeros(1,d);
        slave_i_corner = i_corner;
        counter = 0;
        slave = patch;

        % track counter clockwise around the given vertice for all patches connected to this vertice
        while slave ~= patch || counter == 0
            counter = counter + 1;
            if ~any(isnan(cornerData(slave).v_n)) % skip the cases of "non-Lipschitz vertices"
                v_n = v_n + cornerData(slave).v_n(slave_i_corner,:);
            end
            slave_prev = slave;
            midx = corner2midx(slave_i_corner);
            slave_temp = topologyMap{slave}(midx).slave;
            if isempty(slave_temp)
                idx = mod(find(corner2midx == midx)-2,4)+1;
                midx = corner2midx(idx);
                slave = topologyMap{slave}(midx).slave;
            else
                slave = slave_temp;
            end
            sidx = topologyMap{slave_prev}(midx).sidx;

            slave_i_corner = sidx2corner(sidx);
        end
        cornerData(patch).avg_v_n(i_corner,:) = v_n/counter;
    end
    cornerData(patch).avg_v_n = cornerData(patch).avg_v_n./norm2(cornerData(patch).avg_v_n);
end
% for patch = 1:noBdryPatches
%     quiver3(cornerData(patch).X(:,1), ...
%             cornerData(patch).X(:,2), ...
%             cornerData(patch).X(:,3), ...
%             cornerData(patch).avg_v_n(:,1), ...
%             cornerData(patch).avg_v_n(:,2), ...
%             cornerData(patch).avg_v_n(:,3),'LineWidth',1,'AutoScale','off','DisplayName',['Patch ' num2str(patch)])
% end
%% Compute angle between all pairs of connected patches
midx2corners = [1,4; 2,3; 1,2; 4,3];
for patch = indices
    for midx = 1:2*d_p
        if isnan(angles(patch,midx))
            slave = topologyMap{patch}(midx).slave;
            if isempty(slave)
                continue
            end
            sidx = topologyMap{patch}(midx).sidx;
            orient = topologyMap{patch}(midx).orient;
            sgn_m = (-1)^midx;
            i_corners_m = midx2corners(midx,:);
            switch midx
                case {1,2}
                    v_tm = sgn_m*cornerData(patch).v_t1(i_corners_m,:);
                case {3,4}
                    v_tm = sgn_m*cornerData(patch).v_t2(i_corners_m,:);
            end
            v_nm = cornerData(patch).v_n(i_corners_m,:);
            i_corners_s = midx2corners(sidx,:);
            sgn_s = (-1)^sidx;
            switch sidx
                case {1,2}
                    v_ts = sgn_s*cornerData(slave).v_t1(i_corners_s,:);
                case {3,4}
                    v_ts = sgn_s*cornerData(slave).v_t2(i_corners_s,:);
            end
            v_ns = cornerData(slave).v_n(i_corners_s,:);
            if orient
                v_ts = flipud(v_ts);
                v_ns = flipud(v_ns);
            end
            nanIndices = or(any(isnan(v_tm),2), any(isnan(v_ts),2));
            v_tm(nanIndices,:) = [];
            v_ts(nanIndices,:) = [];
            v_nm(nanIndices,:) = [];
            v_ns(nanIndices,:) = [];
            angle = mean(abs(acos(dot(v_nm,v_ns,2))));
            if all(dot((v_ts+v_tm)/2,(v_ns+v_nm)/2,2) < 0)
                angle = pi - angle;
            else % the boundary at these corners are convex
                angle = pi + angle;
            end
            
            angles(patch,midx) = angle;
            angles(slave,sidx) = angle;
        end
    end
end