function S2Vobj = surfaceToVolume(varargin)
% This routine assumes nurbs_surf to be a closed NURBS-surface with all
% patches having normal vectors pointing outwards.

% Set default values
options = struct('t', 1, ...                                % thickness of patch created from a surface patch having all angles larger than "sharpAngle" w.r.t. neighbouring patches
                 'sharpAngle', 120*pi/180,...               % Threshold for a "sharp" angle
                 'Eps', 1e-10, ...                          % Threshold for assuming two physical points in the l2-norm to be identical
                 'explodeNURBSsurface', true, ...           % Explode patches having C0-elements
                 'S2V_algorithm',{{'A1_1','A13_1'}}, ...    % Specify the desired algorithm to be used. Here 'A1_1' mean ("A"lgorithm for when only side patch 1 is known and the first of the available algorithms)
                 'avg_v_n_threshholdAngle', 11.25*pi/180);  % Threshold angle deviation between the normal vectors v_n avg_v_n

S2Vobj = varargin{1};
if nargin > 1
    if numel(varargin) > 2
        newOptions = varargin(2:end);
    else
        newOptions = varargin{2};
    end
    options = updateOptions(options,newOptions);
end
if isfield(S2Vobj,'nurbs_bdry')
    counter = S2Vobj.counter;
else
    nurbs_surf = cleanNURBS(S2Vobj.nurbs_surf,[],1e-6);
    
    if options.explodeNURBSsurface
        nurbs_surf = explodeNURBS(nurbs_surf);
    end
    d_p = nurbs_surf{1}.d_p;
    
    S2Vobj.noBdryPatches = numel(nurbs_surf);
    if isfield(options,'connection')
        connection = options.connection;
    else
        geometry = getTopology(nurbs_surf);
        connection = geometry.topology.connection;
    end
    S2Vobj.topologyMap = createTopologyMap(connection,S2Vobj.noBdryPatches,d_p);
    
    
    S2Vobj.nurbs_bdry = nurbs_surf;
    S2Vobj = computeCornerData(S2Vobj,options);
    
    
    S2Vobj.nurbs_vol = cell(1,4*numel(nurbs_surf));
    counter = 0;
end

counter = counter + 1;


% Attach new volumetric patch on nurbs_bdry
[S2Vobj.nurbs_vol(counter), nurbs_covered,nurbs_newBdry] = addPatch(S2Vobj,options);

if ~isempty(nurbs_newBdry)
    % Remove faces with zero measure
    zeroMeasure = NURBShasZeroMeasure(nurbs_newBdry);
    nurbs_newBdry(zeroMeasure) = [];
end

if ~isempty(nurbs_newBdry)
    noNewBdryPatches = numel(nurbs_newBdry);
    newBdryPatches = (numel(S2Vobj.nurbs_bdry)+1):(numel(S2Vobj.nurbs_bdry)+noNewBdryPatches);

    % Add new patches to the global set of surface patches, nurbs_bdry
    S2Vobj.nurbs_bdry = [S2Vobj.nurbs_bdry, nurbs_newBdry];

    % Find the surrounding patches to the patches that has been covered
    surroundingPatches = 4*numel(nurbs_covered);
    counter2 = 1;
    for i = 1:numel(nurbs_covered)
        patch_i = nurbs_covered(i);
        for midx_i = 1:numel(S2Vobj.topologyMap{patch_i})
            if ~isempty(S2Vobj.topologyMap{patch_i}(midx_i).slave)
                surroundingPatches(counter2) = S2Vobj.topologyMap{patch_i}(midx_i).slave;
                counter2 = counter2 + 1;
            end
        end
    end
    surroundingPatches(counter2:end) = [];
    surroundingPatches = setdiff(surroundingPatches,nurbs_covered); % Remove the patches that is replaced


    % Update topologyMap
    alteredPatches = [surroundingPatches,newBdryPatches];
    S2Vobj.alteredPatches = alteredPatches;
    geometry = getTopology(S2Vobj.nurbs_bdry(alteredPatches));
    connectionSub = geometry.topology.connection;

    for patch_i = newBdryPatches
        S2Vobj.topologyMap{patch_i} = struct('slave', cell(1, 4), ...
                                             'sidx', cell(1, 4), ...
                                             'orient', cell(1, 4));
    end

    for i = 1:numel(connectionSub)
        patch_i = alteredPatches(connectionSub{i}.Attributes.master);
        slave_i = alteredPatches(connectionSub{i}.Attributes.slave);
        midx_i = connectionSub{i}.Attributes.midx;
        sidx_i = connectionSub{i}.Attributes.sidx;
        orient = connectionSub{i}.Attributes.orient;

        S2Vobj.topologyMap{patch_i}(midx_i).slave = slave_i;
        S2Vobj.topologyMap{patch_i}(midx_i).sidx = sidx_i;
        S2Vobj.topologyMap{patch_i}(midx_i).orient = orient;
    end
    for i = nurbs_covered
        for j = 1:4
            S2Vobj.topologyMap{i}(j).slave = [];
            S2Vobj.topologyMap{i}(j).sidx = [];
            S2Vobj.topologyMap{i}(j).orient = [];
        end
    end

%     try % to update angles
        S2Vobj = computeCornerData(S2Vobj,options,newBdryPatches);
%     catch ME
%         warning(ME.message)
%         S2Vobj.nurbs_vol(counter) = [];
%         counter = counter - 1;
%     end
end
S2Vobj.angles(nurbs_covered,:) = NaN;
S2Vobj.S2Vcompleted = ~(counter == 1 || any(S2Vobj.angles(:) < pi) || any(~isnan(S2Vobj.angles(1:S2Vobj.noBdryPatches,:)),'all'));

if S2Vobj.S2Vcompleted
    S2Vobj.nurbs_vol(counter+1:end) = [];
end

% Collect changes into the "Surface to volume" object
S2Vobj.counter = counter;


function [nurbs, nurbs_covered,nurbs_newBdry] = addPatch(S2Vobj,options,patch,midx,maxNoSharpAngles)
%     [patch,maxNoSharpAngles,minSum] = findNextPatch(original_angles,sharpAngle); % Check for sharp angles in at the initial nurbs_surf first
%     midx = find(original_angles(patch,:) < sharpAngle);
%     if isinf(minSum)
%         [patch,maxNoSharpAngles] = findNextPatch(angles,sharpAngle);
%         midx = find(angles(patch,:) < sharpAngle);
%     end
nurbs_bdry = S2Vobj.nurbs_bdry;
topologyMap = S2Vobj.topologyMap;
cornerData = S2Vobj.cornerData;
edgeData = S2Vobj.edgeData;
angles = S2Vobj.angles;


d = 3; % dimension
t = options.t;
% original_angles = angles;
sharpAngle = options.sharpAngle;
acuteAngleAdjustment = 1.2;
if nargin < 3
    [patch,maxNoSharpAngles] = findNextPatch(angles,sharpAngle,S2Vobj.noBdryPatches);
    midx = find(angles(patch,:) < sharpAngle);
end
faces = cell(1,6);
switch maxNoSharpAngles
    case 0
        faces(1) = nurbs_bdry(patch);
        switch options.S2V_algorithm{maxNoSharpAngles+1}
            case 'A1_1'
                faces(2) = faceFromNormals(patch,S2Vobj,options,t);
            case 'A1_2'
                X = zeros(4,3);
                for i_corner = 1:4
                    X(i_corner,:) = cornerData(1).X(i_corner,:) + t*cornerData(1).v_n(i_corner,:);
                end
                faces(2) = getPrismData('X',X([1,2,4,3],:),'d_p',2);
                faces(1:2) = homogenizeNURBSparametrization(faces(1:2));
        end
        nurbs = loftNURBS({faces(1),faces(2)},1,1);
        nurbs_covered = patch;
        nurbs_newBdry = subNURBS(nurbs,'at',[0,1;1,1;1,1],'outwardPointingNormals',true); 
    case 1 % Loft path along slave
        slave = topologyMap{patch}(midx).slave;
        sidx = topologyMap{patch}(midx).sidx;

        % Set the patch with oposite corner having the least amount of deviation between the vectors v_n and avg_v_n as master patch
        v_n_1 = edgeData(patch).v_n{midx};
        avg_v_n_i = edgeData(patch).avg_v_n{midx};
        angle_v_n_diff_i = mean(acos(abs(dot(v_n_1,avg_v_n_i,2))));

        v_n_1 = edgeData(slave).v_n{midx};
        avg_v_n_i = edgeData(slave).avg_v_n{midx};
        angle_v_n_diff_i_s = mean(acos(abs(dot(v_n_1,avg_v_n_i,2))));
        if angle_v_n_diff_i_s < angle_v_n_diff_i
            patch = slave;
            midx = sidx;
        end
        slave = topologyMap{patch}(midx).slave;
        sidx = topologyMap{patch}(midx).sidx;
        
        % Fix orientation to the standard setup
        faces(1) = nurbs_bdry(patch);
        faces(3) = nurbs_bdry(slave);
        faces{1} = orientNURBS(faces{1}, idx1To_idx2_orient(midx,1)); % Make "midx = 1"
        faces{3}  = orientNURBS(faces{3},  idx1To_idx2_orient(sidx,3));   % Make "sidx = 3"
        knots = [faces{3}.knots(end), faces{1}.knots];
        degree = [faces{3}.degree(end), faces{1}.degree];
        number = [faces{3}.number(end), faces{1}.number];
        master_coeffs = faces{1}.coeffs;
        master_coeffs = reshape(master_coeffs,[d+1,1,faces{1}.number]);
        slave1_coeffs = permute(faces{3}.coeffs,[1,3,2]);
        slave1_coeffs = reshape(slave1_coeffs,[d+1,faces{3}.number(2),1,faces{3}.number(1)]);
        
        switch options.S2V_algorithm{maxNoSharpAngles+1}
            case 'A13_1'
                try
%                     len = vecnorm(slave1_coeffs(1:3,end,1,:) - slave1_coeffs(1:3,1,1,:),2,1);
%                     face_normals = faceFromNormals(patch,S2Vobj,options);
%                     face_normals = orientNURBS(face_normals{1}, idx1To_idx2_orient(midx,1)); % Make "midx = 1"
% 
%                     g = reshape(aveknt(knots{1}, degree(1)+1),1,[]);
%                     coeffs = zeros([d+1,number]);
%                     coeffs(4,:,:,:) = repmat(master_coeffs(4,1,:,:),1,number(1));
%                     coeffs(1:3,:,:,:) = master_coeffs(1:3,:,:,:) + len.*reshape(face_normals.coeffs(1:3,:,:),[d,1,number(2:3)]).*g;
%                     coeffs(:,:,1,:) = slave1_coeffs(:,:,1,:);
%                     nurbs = createNURBSobject(coeffs,knots);
                    face_normals1 = faceFromNormals(patch,S2Vobj,options);
                    face_normals1 = orientNURBS(face_normals1{1}, idx1To_idx2_orient(midx,1)); % Make "midx = 1"
                    face_normals2 = faceFromNormals(slave,S2Vobj,options);
                    face_normals2 = orientNURBS(face_normals2{1}, idx1To_idx2_orient(sidx,3)); % Make "sidx = 3"
                    
                    P1 = master_coeffs(1:3,1,:,:);
                    P3 = slave1_coeffs(1:3,:,1,:);

                    Pdiff1 = repmat(P3 - slave1_coeffs(1:3,1,1,:), 1,1,number(2));
                    Pdiff2 = repmat(P1 - master_coeffs(1:3,1,1,:), 1,number(1));
                    len1 = vecnorm(Pdiff1,2,1);
                    len2 = vecnorm(Pdiff2,2,1);
                    w = len1./(len1+len2);
                    w(:,1,1,:) = 1;

                    % Construct rotation matrix for the master face
                    d1 = Pdiff1(:,2,:,:);
                    d1_hat = d1./vecnorm(d1,2,1);
                    n_tilde1 = reshape(face_normals1.coeffs(1:3,:,:),[d,1,number(2:3)]);
                    u = cross(d1_hat, n_tilde1, 1);
                    theta = abs(acos(dot(d1_hat, n_tilde1,1)));
                    u1 = u(1,:,:,:);
                    u2 = u(2,:,:,:);
                    u3 = u(3,:,:,:);
                    sinTheta = sin(theta);
                    cosTheta = cos(theta);
                    R1row1 = repmat(cat(1, cosTheta + u1.^2.*(1-cosTheta), ...
                                    u1.*u2.*(1-cosTheta) - u3.*sinTheta, ...
                                    u1.*u3.*(1-cosTheta) + u2.*sinTheta), 1,number(1));
                    R1row2 = repmat(cat(1, u2.*u1.*(1-cosTheta) + u3.*sinTheta, ...
                                    cosTheta + u2.^2.*(1-cosTheta), ...
                                    u2.*u3.*(1-cosTheta) - u1.*sinTheta), 1,number(1));
                    R1row3 = repmat(cat(1, u3.*u1.*(1-cosTheta) - u2.*sinTheta, ...
                                    u3.*u2.*(1-cosTheta) + u1.*sinTheta, ...
                                    cosTheta + u3.^2.*(1-cosTheta)), 1,number(1));

                    % Construct rotation matrix for the slave face
                    d2 = Pdiff2(:,:,2,:);
                    d2_hat = d2./vecnorm(d2,2,1);
                    n_tilde2 = reshape(face_normals2.coeffs(1:3,:,:),[d,number(1),1,number(3)]);
                    u = cross(d2_hat, n_tilde2,1);
                    theta = abs(acos(dot(d2_hat, n_tilde2,1)));
                    u1 = u(1,:,:,:);
                    u2 = u(2,:,:,:);
                    u3 = u(3,:,:,:);
                    sinTheta = sin(theta);
                    cosTheta = cos(theta);
                    R2row1 = repmat(cat(1, cosTheta + u1.^2.*(1-cosTheta), ...
                                    u1.*u2.*(1-cosTheta) - u3.*sinTheta, ...
                                    u1.*u3.*(1-cosTheta) + u2.*sinTheta), 1,1,number(2));
                    R2row2 = repmat(cat(1, u2.*u1.*(1-cosTheta) + u3.*sinTheta, ...
                                    cosTheta + u2.^2.*(1-cosTheta), ...
                                    u2.*u3.*(1-cosTheta) - u1.*sinTheta), 1,1,number(2));
                    R2row3 = repmat(cat(1, u3.*u1.*(1-cosTheta) - u2.*sinTheta, ...
                                    u3.*u2.*(1-cosTheta) + u1.*sinTheta, ...
                                    cosTheta + u3.^2.*(1-cosTheta)), 1,1,number(2));

                    coeffs = zeros([d+1,number]);
                    if any(isnan(d1_hat(:))) && any(isnan(d2_hat(:)))
                        error('Cannot perform algorithm on patches both having collapsed faces')
                    elseif any(isnan(d1_hat(:)))
                        coeffs(1:3,:,:,:) =           P3 + cat(1, dot(R2row1,Pdiff2,1), dot(R2row2,Pdiff2,1), dot(R2row3,Pdiff2,1));
                    elseif any(isnan(d2_hat(:)))
                        coeffs(1:3,:,:,:) =           P1 + cat(1, dot(R1row1,Pdiff1,1), dot(R1row2,Pdiff1,1), dot(R1row3,Pdiff1,1));
                    else
                        coeffs(1:3,:,:,:) =       w.*(P1 + cat(1, dot(R1row1,Pdiff1,1), dot(R1row2,Pdiff1,1), dot(R1row3,Pdiff1,1))) ...
                                            + (1-w).*(P3 + cat(1, dot(R2row1,Pdiff2,1), dot(R2row2,Pdiff2,1), dot(R2row3,Pdiff2,1)));
                    end
                    coeffs(1:3,1,:,:) = P1;
                    coeffs(1:3,:,1,:) = P3;
                    coeffs(4,:,:,:) = repmat(master_coeffs(4,1,:,:),1,number(1));

                    nurbs = createNURBSobject(coeffs,knots);
%                     
                    if checkOrientation(nurbs, 10)
                        error('Self intersection!')
                    end
                catch ME
                    warning(ME.message)
                    coeffs = master_coeffs + slave1_coeffs - slave1_coeffs(:,1,1,:);
                    nurbs = createNURBSobject(coeffs,knots);
                end
            case 'A13_2'
                coeffs = master_coeffs + slave1_coeffs - slave1_coeffs(:,1,1,:);
                nurbs = createNURBSobject(coeffs,knots);
        end
        nurbs_covered = [patch,slave];
        nurbs_newBdry = subNURBS(nurbs,'at',[0,1;0,1;1,1],'outwardPointingNormals',true);
    case 2
        midx_regular = singularFromRegular(topologyMap,patch,setdiff(1:4,midx));
        if ~isempty(midx_regular)
            [nurbs, nurbs_covered,nurbs_newBdry] = addPatch(S2Vobj,options,patch,sort([midx,midx_regular]),maxNoSharpAngles+1);
            return
        end
        if all(ceil(midx/2) == [1,2]) % Fill a patch in a corner
            % Set the patch with oposite corner having the least amount of deviation between the vectors v_n and avg_v_n as master patch
            angle_v_n_diff = zeros(1,3);
            midx = correct_midx_order(midx);
            idx_oposite = get_oposite_corner_index(intermediateCorner(midx));
            
            v_n_1 = cornerData(patch).v_n(idx_oposite,:);
            avg_v_n_i = cornerData(patch).avg_v_n(idx_oposite,:);
            angle_v_n_diff(1) = acos(abs(dot(v_n_1,avg_v_n_i)));

            slave3 = topologyMap{patch}(midx(1)).slave;
            sidx3 = topologyMap{patch}(midx(1)).sidx;
            sidx3_m = sort([sidx3,midxToLeftmidx(sidx3)]);
            sidx3_m = correct_midx_order(sidx3_m);
            if any(isnan(angles(slave3,sidx3_m))) % slave1 and slave2 are "connected" through a singularity and patch should remain the master patch
                angle_v_n_diff(2) = Inf;
            else
                idx_oposite = get_oposite_corner_index(intermediateCorner(sidx3_m));
                v_n_1 = cornerData(slave3).v_n(idx_oposite,:);
                avg_v_n_i = cornerData(slave3).avg_v_n(idx_oposite,:);
                angle_v_n_diff(2) = acos(abs(dot(v_n_1,avg_v_n_i)));
            end

            slave5 = topologyMap{patch}(midx(2)).slave;
            sidx5 = topologyMap{patch}(midx(2)).sidx;
            sidx5_m = sort([sidx5,midxToRightmidx(sidx5)]);
            sidx5_m = correct_midx_order(sidx5_m);
            if any(isnan(angles(slave5,sidx5_m))) % slave1 and slave2 are "connected" through a singularity and patch should remain the master patch
                angle_v_n_diff(3) = Inf;
            else
                idx_oposite = get_oposite_corner_index(intermediateCorner(sidx5_m));
                v_n_1 = cornerData(slave5).v_n(idx_oposite,:);
                avg_v_n_i = cornerData(slave5).avg_v_n(idx_oposite,:);
                angle_v_n_diff(3) = acos(abs(dot(v_n_1,avg_v_n_i)));
            end

            [~,I] = min(angle_v_n_diff);
            switch I
                case 2
                    patch = slave3;
                    midx = sidx3_m;
                case 3
                    patch = slave5;
                    midx = sidx5_m;
            end

        
            slave3 = topologyMap{patch}(midx(1)).slave;
            sidx3 = topologyMap{patch}(midx(1)).sidx;

            slave5 = topologyMap{patch}(midx(2)).slave;
            sidx5 = topologyMap{patch}(midx(2)).sidx;
            if isempty(slave3)
                slave3 = topologyMap{slave5}(midxToRightmidx(sidx5)).slave;
                sidx3 = midxToRightmidx(topologyMap{slave3}(midxToRightmidx(sidx5)).sidx);
            end

            if isempty(slave5)
                slave5 = topologyMap{slave3}(midxToLeftmidx(sidx3)).slave;
                sidx5 = midxToLeftmidx(topologyMap{slave3}(midxToLeftmidx(sidx3)).sidx);
            end

            % Fix orientation to the standard setup
            faces(1) = nurbs_bdry(patch);
            faces(3) = nurbs_bdry(slave3);
            faces(5) = nurbs_bdry(slave5);

            faces{1} = orientNURBS(faces{1}, idx1To_idx2_orient(midx(1),1)); % Make "midx = 1"
            faces{3} = orientNURBS(faces{3}, idx1To_idx2_orient(sidx3,3));   % Make "sidx = 3"
            faces{5} = orientNURBS(faces{5}, idx1To_idx2_orient(sidx5,1));   % Make "sidx = 1"


            knots = [faces{3}.knots(end), faces{1}.knots];
            degree = [faces{3}.degree(end), faces{1}.degree];
            number = [faces{3}.number(end), faces{1}.number];
            master_coeffs = faces{1}.coeffs;
            master_coeffs = reshape(master_coeffs,[d+1,1,faces{1}.number]);
            slave1_coeffs = permute(faces{3}.coeffs,[1,3,2]);
            slave1_coeffs = reshape(slave1_coeffs,[d+1,faces{3}.number(2),1,faces{3}.number(1)]);
            slave2_coeffs = faces{5}.coeffs;

            try
%                 midx_opposite = midx-(-1).^midx;
%                 cornerIndices = intermediateCorner(midx);
%                 oppositeCornerIndices = intermediateCorner(midx_opposite);
% 
%                 normals = cornerData(patch).avg_v_n(oppositeCornerIndices,:).';
%                 v_n = cornerData(patch).v_n(oppositeCornerIndices,:).';
%                 indices = acos(abs(dot(normals,v_n,1))) > avg_v_n_threshholdAngle;
%                 normals(:, indices) = v_n(:,indices);
%                 acuteConnection = (angles(patch,midx) + angles(patch,midx_opposite)/2) < pi;
%                 if true
%                     if any(acuteConnection)
%                         error('acute angles found')
%                     end
%                 else
%                     if all(acuteConnection) % Make opposite sides "parallel"
%                         normals = cornerData(patch).v_t2(cornerIndices,:).';
%                     else
%                         for i = 1:2
%                             sgn = (-1)^midx(i);
%                             if acuteConnection(i)
%                                 switch midx(i)
%                                     case {1,2}
%                                         v_t = cornerData(patch).v_xi(oppositeCornerIndices,:);
%                                     case {3,4}
%                                         v_t = cornerData(patch).v_eta(oppositeCornerIndices,:);
%                                 end
%                                 v_n = cornerData(patch).v_n(oppositeCornerIndices,:);
% %                                 phi = acos(dot(v_t,v_n));
%                                 theta = angles(patch,midx(i))*acuteAngleAdjustment;
% %                                 cornerNormal = ([v_t; v_n]\[cos(theta); cos(theta-phi)]);
% %                                 cornerNormal = cornerNormal./vecnorm(cornerNormal,1);
%         
%                                 normals = (-sgn*v_t*cos(theta) + v_n*sin(theta)).';
%                             end
%                         end
%                     end
%                 end
% 
%                 len2 = norm(slave1_coeffs(1:3,end,1,end) - slave1_coeffs(1:3,1,1,end));
%                 g = reshape(aveknt(knots{1}, degree(1)+1),1,[]);
%                 coeffs2 = zeros(d+1,number(1));
%                 coeffs2(1:3,:) = master_coeffs(1:3,1,end,end) + len2*normals.*g;
%                 coeffs2(4,:) = [master_coeffs(4,1,end,end), slave1_coeffs(4,2:end,1,end)]; % Use the weights based on the weights on opposite side
% 
%                 edges = cell(1,4);
%                 edges(1) = subNURBS(faces(5),'at',[0,0;0,1]);
%                 edges(2) = createNURBSobject(coeffs2,knots{1});
%                 edges(3) = subNURBS(faces(1),'at',[0,1;0,0]);
%                 g = reshape(aveknt(knots{3}, degree(3)+1),1,[]);
%                 coeffs4 = edges{1}.coeffs(:,end).*(1-g) + coeffs2(:,end).*g;
%                 edges(4) = createNURBSobject(coeffs4,knots{3});
%                 faces(4) = GordonHall(edges);
% 
%                 edges(1) = subNURBS(faces(3),'at',[0,0;0,1]);
%                 edges(2) = subNURBS(faces(4),'at',[0,0;0,1]);
%                 edges(3) = subNURBS(faces(5),'at',[0,1;0,0]);
%                 g = reshape(aveknt(knots{2}, degree(2)+1),1,[]);
%                 coeffs4 = edges{1}.coeffs(:,end).*(1-g) + edges{2}.coeffs(:,end).*g;
%                 edges(4) = createNURBSobject(coeffs4,knots{2});
%                 faces(2) = GordonHall(edges);
% 
%                 edges(1) = subNURBS(faces(1),'at',[0,0;0,1]);
%                 edges(2) = subNURBS(faces(2),'at',[0,0;0,1]);
%                 edges(3) = subNURBS(faces(3),'at',[0,1;0,0]);
%                 edges(4) = subNURBS(faces(4),'at',[0,1;0,0]);
%                 faces(6) = GordonHall(edges);
%                 
%                 nurbs = GordonHall(faces);
%                 faces_old = faces;
%                 opposite_patchID = [patch,slave3,slave5];
%                 counter = 1;
%                 for i_face = [1,3,5]
%                     if checkOrientation(nurbs, 10)
%                         faces = faces_old;
%                         face = faceFromNormals(opposite_patchID(counter),S2Vobj,options,len2);
%                         face = orientNURBS(face{1}, idx1To_idx2_orient(midx,1)); % Make "midx = 1"
%                         faces{i_face+1}.coeffs(:,2:end,2:end) = face{1}.coeffs(:,2:end,2:end);
%                         switch i_face
%                             case 1
%                                 faces{4}.coeffs(:,2:end,end) = face{1}.coeffs(:,end,2:end);
%                                 faces{6}.coeffs(:,end,2:end) = face{1}.coeffs(:,2:end,end);
%                             case 3
%                                 faces{6}.coeffs(:,2:end,end) = face{1}.coeffs(:,end,2:end);
%                                 faces{2}.coeffs(:,end,2:end) = face{1}.coeffs(:,2:end,end);
%                             case 5
%                                 faces{2}.coeffs(:,2:end,end) = face{1}.coeffs(:,end,2:end);
%                                 faces{4}.coeffs(:,end,2:end) = face{1}.coeffs(:,2:end,end);
%                         end
%                         nurbs = GordonHall(faces);
%                     else
%                         break
%                     end
%                     counter = counter + 1;
%                 end
                if 1
                    len1 = repmat(vecnorm(slave1_coeffs(1:3,end,1,:) - slave1_coeffs(1:3,1,1,:),2,1), [1,1,number(2)]);
                    len2 = repmat(vecnorm(slave2_coeffs(1:3,end,:,1) - slave2_coeffs(1:3,1,:,1),2,1), [1,1,1,number(3)]);
                    len = mean(cat(1,len1,len2),1);
                    face_normals = faceFromNormals(patch,S2Vobj,options);
                    face_normals = orientNURBS(face_normals{1}, idx1To_idx2_orient(midx,1)); % Make "midx = 1"
    
                    g = reshape(aveknt(knots{1}, degree(1)+1),1,[]);
                    coeffs = zeros([d+1,number]);
                    coeffs(4,:,:,:) = repmat(master_coeffs(4,1,:,:),1,number(1));
                    coeffs(1:3,:,:,:) = master_coeffs(1:3,:,:,:) + len.*reshape(face_normals.coeffs(1:3,:,:),[d,1,number(2:3)]).*g;
                    coeffs(:,:,1,:) = slave1_coeffs(:,:,1,:);
                    coeffs(:,:,:,1) = slave2_coeffs(:,:,:,1);
                    nurbs = createNURBSobject(coeffs,knots);
                else
                    len1 = repmat(vecnorm(slave1_coeffs(1:3,end,1,:) - slave1_coeffs(1:3,1,1,:),2,1), [1,1,number(2)]);
                    len2 = repmat(vecnorm(slave2_coeffs(1:3,end,:,1) - slave2_coeffs(1:3,1,:,1),2,1), [1,1,1,number(3)]);
                    len = mean(cat(1,len1,len2),1);
                    face_normals = faceFromNormals(patch,S2Vobj,options);
                    face_normals = orientNURBS(face_normals{1}, idx1To_idx2_orient(midx,1)); % Make "midx = 1"
    
                    g = reshape(aveknt(knots{1}, degree(1)+1),1,[]);
                    coeffs = zeros([d+1,number]);
                    coeffs(4,:,:,:) = repmat(master_coeffs(4,1,:,:),1,number(1));
                    coeffs(1:3,:,:,:) = master_coeffs(1:3,:,:,:) + len.*reshape(face_normals.coeffs(1:3,:,:),[d,1,number(2:3)]).*g;
                    coeffs(:,:,1,:) = slave1_coeffs(:,:,1,:);
                    coeffs(:,:,:,1) = slave2_coeffs(:,:,:,1);
                    nurbs = createNURBSobject(coeffs,knots);

%                     face_normals1 = faceFromNormals(patch,S2Vobj,options);
%                     face_normals1 = orientNURBS(face_normals1{1}, idx1To_idx2_orient(midx,1)); % Make "midx = 1"
%                     face_normals2 = faceFromNormals(slave3,S2Vobj,options);
%                     face_normals2 = orientNURBS(face_normals2{1}, idx1To_idx2_orient(sidx3,3)); % Make "sidx = 3"
%                     face_normals3 = faceFromNormals(slave5,S2Vobj,options);
%                     face_normals3 = orientNURBS(face_normals3{1}, idx1To_idx2_orient(sidx5,1)); % Make "sidx = 1"
% 
%                     g1 = reshape(aveknt(knots{1}, degree(1)+1),1,[]);
%                     g2 = reshape(aveknt(knots{2}, degree(1)+1),1,1,[]);
%                     g3 = reshape(aveknt(knots{3}, degree(1)+1),1,1,1,[]);
%                     coeffs =   (master_coeffs(:,1,:,:) + len1.*face_normals1).*g1/3 ...
%                              + (slave1_coeffs(:,:,1,:) + len2.*face_normals2).*g2/3 ...
%                              + (slave2_coeffs(:,:,:,1) + len3.*face_normals3).*g3/3 ...
%                              + master_coeffs(:,1,1,1).*g1.*g2.*g3;
% 
%                     coeffs = zeros([d+1,number]);
%                     edgeCoeff3
%                     avg1 = mean(cat(edgeCoeff3,edgeCoeff5));
%                     ilen = 1./len1 + 1./len2 + 1./len3;
%                     coeffs(1:3,end,end,end) =   (master_coeffs(:,1,end,end) + len1.*face_normals1).*1./len1./ilen ...
%                                               + (slave1_coeffs(:,end,1,end) + len2.*face_normals2).*1./len2./ilen ...
%                                               + (slave2_coeffs(:,end,end,1) + len3.*face_normals3).*1./len3./ilen;

                end
                if checkOrientation(nurbs, 10)
                    error('Self intersection!')
                end
            catch ME
                warning(ME.message)
                noStrategies = 2;
                nurbs = cell(1,noStrategies);
                mean_J_r = zeros(1,noStrategies);
                xi = linspace(0,1,10);
                [XI,ETA,ZETA] = ndgrid(xi);
                for type = 1:noStrategies
                    switch type
                        case 1
                            coeffs = master_coeffs + slave1_coeffs + slave2_coeffs - master_coeffs(:,1,1,:) - slave2_coeffs(:,1,:,1) - slave1_coeffs(:,:,1,1) + master_coeffs(:,1,1,1);
                        case 2
                            slave3_coeffs = master_coeffs + slave2_coeffs - slave2_coeffs(:,1,:,1);
                            g = reshape(aveknt(faces{1}.knots{1},faces{1}.degree(1)+1),1,1,[]);
                    
                            coeffs = (master_coeffs + slave1_coeffs + slave2_coeffs - master_coeffs(:,1,1,:) - slave2_coeffs(:,1,:,1) - slave1_coeffs(:,:,1,1) + master_coeffs(:,1,1,1)).*(1-g) ...
                                   + (master_coeffs + slave3_coeffs(:,:,end,:) + slave2_coeffs - master_coeffs(:,1,end,:) - slave2_coeffs(:,1,:,1) - slave3_coeffs(:,:,end,1) + master_coeffs(:,1,end,1)).*g;
                    end
                    nurbs(type) = createNURBSobject(coeffs,knots);
                    mean_J_r(type) = mean(meanRatioJacobian(nurbs{type},[XI(:),ETA(:),ZETA(:)]));                    
                end
                [~,I] = max(mean_J_r);
                nurbs = nurbs(I);
            end
            nurbs_covered = [patch,slave3,slave5];
            nurbs_newBdry = subNURBS(nurbs,'at',[0,1;0,1;0,1],'outwardPointingNormals',true);
        else
            slave3 = topologyMap{patch}(midx(1)).slave;
            sidx3 = topologyMap{patch}(midx(1)).sidx;
            slave4 = topologyMap{patch}(midx(2)).slave;
            sidx4 = topologyMap{patch}(midx(2)).sidx;

            faces(1) = nurbs_bdry(patch);
            faces(3) = nurbs_bdry(slave3);
            faces(4) = nurbs_bdry(slave4);

            % Fix orientation to the standard setup
            faces{1} = orientNURBS(faces{1}, idx1To_idx2_orient(midx(1),1)); % Make "midx = 1"
            faces{3} = orientNURBS(faces{3}, idx1To_idx2_orient(sidx3,3));   % Make "sidx = 3"
            faces{4} = orientNURBS(faces{4}, idx1To_idx2_orient(sidx4,3,true));   % Make "sidx = 1"


            knots = [faces{3}.knots(end), faces{1}.knots];
            master_coeffs = faces{1}.coeffs;
            master_coeffs = reshape(master_coeffs,[d+1,1,faces{1}.number]);
            slave1_coeffs = permute(faces{3}.coeffs,[1,3,2]);
            slave1_coeffs = reshape(slave1_coeffs,[d+1,faces{3}.number(2),1,faces{3}.number(1)]);
            slave2_coeffs = permute(faces{4}.coeffs,[1,3,2]);
            slave2_coeffs = reshape(slave2_coeffs,[d+1,faces{4}.number(2),1,faces{4}.number(1)]);
            g = reshape(aveknt(faces{1}.knots{1},faces{1}.degree(1)+1),1,1,[]);
            coeffs = (master_coeffs + slave1_coeffs - slave1_coeffs(:,1,1,:)).*(1-g) + (master_coeffs + slave2_coeffs - slave2_coeffs(:,1,1,:)).*g;

            nurbs = createNURBSobject(coeffs,knots);

            nurbs_covered = [patch,slave3,slave4];
            nurbs_newBdry = subNURBS(nurbs,'at',[0,1;0,0;1,1],'outwardPointingNormals',true);
        end
    case 3
        midx_regular = singularFromRegular(topologyMap,patch,setdiff(1:4,midx));
        if ~isempty(midx_regular)
            [nurbs, nurbs_covered,nurbs_newBdry] = addPatch(S2Vobj,options,patch,sort([midx,midx_regular]),maxNoSharpAngles+1);
            return
        end
        midx = correct_midx_order(midx);
        slave3 = topologyMap{patch}(midx(1)).slave;
        sidx3 = topologyMap{patch}(midx(1)).sidx;
        slave5 = topologyMap{patch}(midx(2)).slave;
        sidx2 = topologyMap{patch}(midx(2)).sidx;
        if isempty(slave5)
            slave_temp = topologyMap{patch}(midxToLeftmidx(midx(2))).slave;
            sidx_temp = topologyMap{patch}(midxToLeftmidx(midx(2))).sidx;
            slave5 = topologyMap{slave_temp}(midxToLeftmidx(sidx_temp)).slave;
            sidx2 = midxToLeftmidx(topologyMap{slave_temp}(midxToLeftmidx(sidx_temp)).sidx);
        end

        slave4 = topologyMap{patch}(midx(3)).slave;
        sidx4 = topologyMap{patch}(midx(3)).sidx;

        faces(1) = nurbs_bdry(patch);
        faces(3) = nurbs_bdry(slave3);
        faces(5) = nurbs_bdry(slave5);
        faces(4) = nurbs_bdry(slave4);

        % Fix orientation to the standard setup
        faces{1} = orientNURBS(faces{1}, idx1To_idx2_orient(midx(1),1)); % Make "midx = 1"
        faces{3} = orientNURBS(faces{3}, idx1To_idx2_orient(sidx3,3));   % Make "sidx = 3"
        faces{5} = orientNURBS(faces{5}, idx1To_idx2_orient(sidx2,1));   % Make "sidx = 1"
        faces{4} = orientNURBS(faces{4}, idx1To_idx2_orient(sidx4,3,true));   % Make "sidx = 1"


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

        nurbs = createNURBSobject(coeffs,knots);

        nurbs_covered = [patch,slave3,slave5,slave4];
        nurbs_newBdry = subNURBS(nurbs,'at',[0,1;0,0;0,1],'outwardPointingNormals',true);
    case 4
        % First check if "maxNoSharpAngles=5" (i.e. all faces are known)
        slavesSlave = zeros(1,4);
        for midx2 = 1:4
            slave = topologyMap{patch}(midx2).slave;
            sidx1 = topologyMap{patch}(midx2).sidx;
            sidx1_opposite = sidx1-(-1)^sidx1;
            slavesSlave(midx2) = topologyMap{slave}(sidx1_opposite).slave;
        end
        slave3 = topologyMap{patch}(midx(1)).slave;
        slave4 = topologyMap{patch}(midx(2)).slave;
        slave5 = topologyMap{patch}(midx(3)).slave;
        slave6 = topologyMap{patch}(midx(4)).slave;
        sidx3 = topologyMap{patch}(midx(1)).sidx;
        sidx4 = topologyMap{patch}(midx(2)).sidx;
        sidx5 = topologyMap{patch}(midx(3)).sidx;
        sidx6 = topologyMap{patch}(midx(4)).sidx;
        faces(1) = nurbs_bdry(patch);

        faces(3) = nurbs_bdry(slave3);
        faces(4) = nurbs_bdry(slave4);
        faces(5) = nurbs_bdry(slave5);
        faces(6) = nurbs_bdry(slave6);

        % Fix orientation to the standard setup
        faces{1} = orientNURBS(faces{1}, idx1To_idx2_orient(midx(1),1));        % Make "midx = 1"
        faces{3} = orientNURBS(faces{3}, idx1To_idx2_orient(sidx3,3));          % Make "sidx = 3"
        faces{5} = orientNURBS(faces{5}, idx1To_idx2_orient(sidx5,1));          % Make "sidx = 1"
        faces{4} = orientNURBS(faces{4}, idx1To_idx2_orient(sidx4,3,true));     % Make "sidx = 3"
        faces{6} = orientNURBS(faces{6}, idx1To_idx2_orient(sidx6,1,true));     % Make "sidx = 1"
        if slavesSlave(1) == slavesSlave(2) && slavesSlave(3) == slavesSlave(4)
            faces(2) = nurbs_bdry(slavesSlave(1));    
            faces{2} = orientNURBS(faces{2}, idx1To_idx2_orient(sidx1,3,true));     % Make "sidx = 3"

            nurbs_covered = [patch,slave3,slave5,slave4,slave6,slavesSlave(1)];
            nurbs_newBdry = [];
        else
            edges = cell(1,4);
            edges(1) = subNURBS(faces{3},'at',[0,0;0,1]);
            edges(2) = subNURBS(faces{4},'at',[0,0;0,1]);
            edges(3) = subNURBS(faces{5},'at',[0,1;0,0]);
            edges(4) = subNURBS(faces{6},'at',[0,1;0,0]);
            faces(2) = GordonHall(edges);
            nurbs_covered = [patch,slave3,slave5,slave4,slave6];
            nurbs_newBdry = faces(2);
        end
        nurbs = GordonHall(faces);

end

function [patch,maxNoSharpAngles,minSum] = findNextPatch(angles,sharpAngle,noBdryPatches)

sharpAngles = angles < sharpAngle;
indices_op = or(and(sharpAngles(:,1),sharpAngles(:,2)),and(sharpAngles(:,3),sharpAngles(:,4)));
indices_covered = all(isnan(angles),2);
indices_nan = and(any(isnan(angles),2),~indices_covered);

noSharpAngles = sum(sharpAngles,2);
noSharpAngles(indices_op) = noSharpAngles(indices_op) + 0.1; % Prioritize cases where oposite faces are known
noSharpAngles(indices_nan) = noSharpAngles(indices_nan) + 0.2; % Prioritize cases containing singularities
noSharpAngles(1:noBdryPatches) = noSharpAngles(1:noBdryPatches) + 0.21; % Prioritize input surface
maxNoSharpAngles = max(noSharpAngles);
candidates = and(noSharpAngles == maxNoSharpAngles,~indices_covered);
avgAngles = Inf(size(sharpAngles,1),1);
for i = 1:size(sharpAngles,1)
    if candidates(i)
        angles_i = angles(i,:);
        angles_i(isnan(angles_i)) = [];
        avgAngles(i) = min(angles_i);
    end
end
maxNoSharpAngles = floor(maxNoSharpAngles);
[minSum,patch] = min(avgAngles);

function S2Vobj = computeCornerData(S2Vobj,options,indices)

nurbs_bdry = S2Vobj.nurbs_bdry;
noBdryPatches = numel(nurbs_bdry);
d_p = nurbs_bdry{1}.d_p;
d = 3; % dimension
if nargin < 3
    indices = 1:noBdryPatches;
end
if ~isfield(S2Vobj,'cornerData')
    S2Vobj.cornerData = struct('X', cell(1, noBdryPatches), ...
                        'v_t1', cell(1, noBdryPatches), ...
                        'v_t2', cell(1, noBdryPatches), ...
                        'v_xi', cell(1, noBdryPatches), ...
                        'v_eta', cell(1, noBdryPatches), ...
                        'cornerAngle', cell(1, noBdryPatches), ...
                        'v_n', cell(1, noBdryPatches), ...
                        'avg_v_n', cell(1, noBdryPatches));
    S2Vobj.edgeData = struct('X', cell(1, noBdryPatches), ...
                        'v_xi', cell(1, noBdryPatches), ...
                        'v_eta', cell(1, noBdryPatches), ...
                        'v_n', cell(1, noBdryPatches), ...
                        'avg_v_n', cell(1, noBdryPatches), ...
                        'collapsed', cell(1, noBdryPatches));
    S2Vobj.angles = NaN(noBdryPatches,2*d_p);
else
    noNewFields = noBdryPatches-numel(S2Vobj.cornerData);
    S2Vobj.cornerData = [S2Vobj.cornerData, struct('X', cell(1, noNewFields), ...
                                     'v_t1', cell(1, noNewFields), ...
                                     'v_t2', cell(1, noNewFields), ...
                                     'v_xi', cell(1, noNewFields), ...
                                     'v_eta', cell(1, noNewFields), ...
                                     'cornerAngle', cell(1, noNewFields), ...
                                     'v_n', cell(1, noNewFields), ...
                                     'avg_v_n', cell(1, noNewFields))];
    S2Vobj.edgeData = [S2Vobj.edgeData, struct('X', cell(1, noNewFields), ...
                                     'v_xi', cell(1, noNewFields), ...
                                     'v_eta', cell(1, noNewFields), ...
                                     'v_n', cell(1, noNewFields), ...
                                     'avg_v_n', cell(1, noNewFields), ...
                                     'collapsed', cell(1, noNewFields))];
    S2Vobj.angles = [S2Vobj.angles; NaN(noNewFields,2*d_p)];
end
Eps = options.Eps;


xi = [0,0;
      1,0;
      1,1;
      0,1]; % Corner indices are oriented counter-clockwise starting at xi,eta=0

%% Compute normal vector at the vertices
for patch = indices
    [X,dXdxi,dXdeta] = evaluateNURBSvec(nurbs_bdry{patch},xi,1);
    v_xi = dXdxi./norm2(dXdxi);
    v_eta = dXdeta./norm2(dXdeta);
    v_t1 = [v_xi(1,:); v_eta(2,:); -v_xi(3,:); -v_eta(4,:)];
    v_t2 = [v_eta(1,:); -v_xi(2,:); -v_eta(3,:); v_xi(4,:)];

    % Handle edges collapsed into a single point
    if norm(X(1,:)-X(2,:)) < Eps
        v_t1(1,:) = v_t1(2,:);
        v_t2(2,:) = v_t2(1,:);
    end
    if norm(X(2,:)-X(3,:)) < Eps
        v_t1(2,:) = v_t1(3,:);
        v_t2(3,:) = v_t2(2,:);
    end
    if norm(X(3,:)-X(4,:)) < Eps
        v_t1(3,:) = v_t1(4,:);
        v_t2(4,:) = v_t2(3,:);
    end
    if norm(X(4,:)-X(1,:)) < Eps
        v_t1(4,:) = v_t1(1,:);
        v_t2(1,:) = v_t2(4,:);
    end

    % Approximate the normal vectors by shifting the evaluation point slightly away from singularity
    % This part could probably be improved
    Eps2 = 1e-5;
    v_n = cross(v_t1, v_t2, 2);
    singularCorners_v_n = find(or(norm2(v_n) < Eps,any(isnan(v_n),2)));
    singularCorners = find(or(norm2(v_n) < Eps, or(norm2(dXdxi) < Eps,norm2(dXdeta) < Eps)));
    if ~isempty(singularCorners)
        xi_eps = [Eps2,Eps2;
                  1-Eps2,Eps2;
                  1-Eps2,1-Eps2;
                  Eps2,1-Eps2]; 
        [~,dXdxi_m2,dXdeta_m2] = evaluateNURBSvec(nurbs_bdry{patch},xi_eps(singularCorners,:),1);
        v_xi(singularCorners,:) = dXdxi_m2./norm2(dXdxi_m2);
        v_eta(singularCorners,:) = dXdeta_m2./norm2(dXdeta_m2);
        if ~isempty(singularCorners_v_n)
            v_n(singularCorners_v_n,:) = cross(v_xi(singularCorners_v_n,:), v_eta(singularCorners_v_n,:), 2);
        end
    end
    cornerAngle = acos(dot(v_t1,v_t2,2));

    v_n = v_n./norm2(v_n);
    S2Vobj.cornerData(patch).X = X;
    S2Vobj.cornerData(patch).v_xi = v_xi;
    S2Vobj.cornerData(patch).v_eta = v_eta;
    S2Vobj.cornerData(patch).v_t1 = v_t1;
    S2Vobj.cornerData(patch).v_t2 = v_t2;
    S2Vobj.cornerData(patch).v_n = v_n;
    S2Vobj.cornerData(patch).cornerAngle = cornerAngle;
    if any(isnan(S2Vobj.cornerData(patch).v_n(:)))
        error('Could not compute normal vector at corner. Case not accounted for')
    end
end

%% Compute normal vector at the edges
for patch = indices
    knots = nurbs_bdry{patch}.knots;
    degree = nurbs_bdry{patch}.degree;
    number = nurbs_bdry{patch}.number;
    S2Vobj.edgeData(patch).collapsed = zeros(1,4);
    S2Vobj.edgeData(patch).X = cell(1,4);
    S2Vobj.edgeData(patch).v_xi = cell(1,4);
    S2Vobj.edgeData(patch).v_eta = cell(1,4);
    S2Vobj.edgeData(patch).v_n = cell(1,4);
    S2Vobj.edgeData(patch).avg_v_n = cell(1,4);
    for midx = 1:4
        switch midx
            case 1
                xi = [zeros(number(2),1), aveknt(knots{2},degree(2)+1).'];
            case 2
                xi = [ones(number(2),1), aveknt(knots{2},degree(2)+1).'];
            case 3
                xi = [aveknt(knots{1},degree(1)+1).', zeros(number(1),1)];
            case 4
                xi = [aveknt(knots{1},degree(1)+1).', ones(number(1),1)];
        end
        [X,dXdxi,dXdeta] = evaluateNURBSvec(nurbs_bdry{patch},xi,1);
        v_xi = dXdxi./norm2(dXdxi);
        v_eta = dXdeta./norm2(dXdeta);
        v_n = cross(v_xi, v_eta, 2);
    
        % Handle edges collapsed into a single point
        if norm(X(1,:)-X(end,:)) < Eps
            S2Vobj.edgeData(patch).collapsed(midx) = true;
            v_n = repmat(S2Vobj.cornerData(patch).v_n(midx2leftCorner(midx),:),size(v_n,1),1);
        else
            S2Vobj.edgeData(patch).collapsed(midx) = false;
            i_corners = midx2corners(midx);
            v_n(1,:) = S2Vobj.cornerData(patch).v_n(i_corners(1),:);
            v_n(end,:) = S2Vobj.cornerData(patch).v_n(i_corners(2),:);
        end
        v_n = v_n./norm2(v_n); % Normalize normal vector

        S2Vobj.edgeData(patch).X{midx} = X;
        S2Vobj.edgeData(patch).v_xi{midx} = v_xi;
        S2Vobj.edgeData(patch).v_eta{midx} = v_eta;
        S2Vobj.edgeData(patch).v_n{midx} = v_n;
        S2Vobj.edgeData(patch).avg_v_n{midx} = NaN(size(v_n)); % Allocation for later
    end
end

%% Compute average normal vector in all vertices
noCorners = 2^d_p;
for patch = indices
    S2Vobj.cornerData(patch).avg_v_n = zeros(4,d);
    for i_corner = 1:noCorners
        avg_v_n = zeros(1,d);
        slave_i_corner = i_corner;
        counter = 0;
        slave = patch;
        
        neighbours = zeros(100,2); % Allocate redundant amount of neighbours
        % track counter clockwise around the given vertice for all patches connected to this vertice
        while slave ~= patch || counter == 0
            counter = counter + 1;
            if counter > 100
                error('Topology is inconsistent')
            end

            if ~any(isnan(S2Vobj.cornerData(slave).v_n))
                cornerAngle = S2Vobj.cornerData(slave).cornerAngle(slave_i_corner);
                avg_v_n = avg_v_n + cornerAngle*S2Vobj.cornerData(slave).v_n(slave_i_corner,:)/(2*pi); % weight the normal vector with the angle of the respetive corner
            end
            slave_prev = slave;
            midx = corner2leftmidx(slave_i_corner);
            slave_temp = S2Vobj.topologyMap{slave}(midx).slave;
            if isempty(slave_temp)
                idx = mod(find(corner2leftmidx(1:4) == midx)-2,4)+1;
                midx = corner2leftmidx(idx);
                slave = S2Vobj.topologyMap{slave}(midx).slave;
            else
                slave = slave_temp;
            end
            sidx = S2Vobj.topologyMap{slave_prev}(midx).sidx;

            slave_i_corner = midx2leftCorner(sidx);
            neighbours(counter,1) = slave;
            neighbours(counter,2) = slave_i_corner;
        end
        neighbours(counter+1,1) = patch;
        neighbours(counter+1,2) = i_corner;
        avg_v_n = avg_v_n./norm2(avg_v_n); % Normalize normal vector
        
        % Update the average normal vector for all patches sharing this corner
        for i = 1:counter+1
            slave = neighbours(i,1);
            slave_i_corner = neighbours(i,2);
            X = S2Vobj.cornerData(slave).X(slave_i_corner,:);
            for j = 1:4
                if norm(X-S2Vobj.cornerData(slave).X(j,:)) < Eps
                    S2Vobj.cornerData(slave).avg_v_n(j,:) = avg_v_n;
%                     sidx = corner2leftmidx(j);
                    sidx = corners2midxes(j);
                    for k = 1:numel(sidx)
                        slave2 = S2Vobj.topologyMap{slave}(sidx(k)).slave;
                        if isempty(slave2)
                            S2Vobj.edgeData(slave).avg_v_n{sidx(k)} = repmat(avg_v_n, size(S2Vobj.edgeData(slave).avg_v_n{sidx(k)},1),1);
                        else
                            if ((sidx(k) == 2 || sidx(k) == 3) && k == 1) || ((sidx(k) == 1 || sidx(k) == 4) && k == 2)
                                S2Vobj.edgeData(slave).avg_v_n{sidx(k)}(end,:) = avg_v_n;
                            else
                                S2Vobj.edgeData(slave).avg_v_n{sidx(k)}(1,:) = avg_v_n;
                            end
%                             if ~isempty(slave2)
%                                 sidx2 = S2Vobj.topologyMap{slave}(sidx(k)).sidx;
%                                 if sidx2 == 1 || sidx2 == 4
%                                     S2Vobj.edgeData(slave2).avg_v_n{sidx2}(end,:) = avg_v_n;
%                                 else
%                                     S2Vobj.edgeData(slave2).avg_v_n{sidx2}(1,:) = avg_v_n;
%                                 end
%                                 if ((sidx2(k) == 1 || sidx2(k) == 4) && k == 1) || ((sidx2(k) == 2 || sidx2(k) == 3) && k == 2)
%                                     S2Vobj.edgeData(slave2).avg_v_n{sidx2(k)}(end,:) = avg_v_n;
%                                 else
%                                     S2Vobj.edgeData(slave2).avg_v_n{sidx2(k)}(1,:) = avg_v_n;
%                                 end
%                             end
                        end
                    end
                end
            end
        end
    end
end

%% Compute remaining average normal vector at the edges (2:end-1) on all edges
for patch = indices
    for midx = 1:4
        slave = S2Vobj.topologyMap{patch}(midx).slave;
        if isempty(slave) % edge midx is collapsed
            S2Vobj.edgeData(patch).avg_v_n{midx} = repmat(S2Vobj.cornerData(patch).avg_v_n(midx2leftCorner(midx),:),...
                                                                        size(S2Vobj.edgeData(patch).v_n{midx},1), 1);
        else
            sidx = S2Vobj.topologyMap{patch}(midx).sidx;
            v_nm = S2Vobj.edgeData(patch).v_n{midx}(2:end-1,:);
            v_ns = S2Vobj.edgeData(slave).v_n{sidx}(2:end-1,:);
            avg_v_n = (v_nm+v_ns)/2;
%             i_corners = midx2corners(midx);
%             avg_v_n(1,:) = S2Vobj.cornerData(patch).avg_v_n(i_corners(1),:);
%             avg_v_n(end,:) = S2Vobj.cornerData(patch).avg_v_n(i_corners(2),:);
            S2Vobj.edgeData(patch).avg_v_n{midx}(2:end-1,:) = avg_v_n./norm2(avg_v_n);
            S2Vobj.edgeData(slave).avg_v_n{sidx}(2:end-1,:) = S2Vobj.edgeData(patch).avg_v_n{midx}(2:end-1,:); % Update the slave patch
        end
    end
end
for patch = 1:indices(end)
    for midx = 1:4
        if any(isnan(S2Vobj.edgeData(patch).avg_v_n{midx}))
            error('NaN encountered - Something is not implemented correctly!')
        end
    end
end


%% Compute angle between all pairs of connected patches
for patch = indices
    for midx = 1:2*d_p
        if isnan(S2Vobj.angles(patch,midx))
            slave = S2Vobj.topologyMap{patch}(midx).slave;
            if isempty(slave)
                continue
            end
            sidx = S2Vobj.topologyMap{patch}(midx).sidx;
            orient = S2Vobj.topologyMap{patch}(midx).orient;
            sgn_m = (-1)^midx;
            if false
                i_corners_m = midx2corners(midx);
                switch midx
                    case {1,2}
                        v_tm = sgn_m*S2Vobj.cornerData(patch).v_xi(i_corners_m,:);
                    case {3,4}
                        v_tm = sgn_m*S2Vobj.cornerData(patch).v_eta(i_corners_m,:);
                end
                v_nm = S2Vobj.cornerData(patch).v_n(i_corners_m,:);
                i_corners_s = midx2corners(sidx);
                sgn_s = (-1)^sidx;
                switch sidx
                    case {1,2}
                        v_ts = sgn_s*S2Vobj.cornerData(slave).v_xi(i_corners_s,:);
                    case {3,4}
                        v_ts = sgn_s*S2Vobj.cornerData(slave).v_eta(i_corners_s,:);
                end
                v_ns = S2Vobj.cornerData(slave).v_n(i_corners_s,:);
                if orient
                    v_ts = flipud(v_ts);
                    v_ns = flipud(v_ns);
                end
            else
                switch midx
                    case {1,2}
                        v_tm = sgn_m*S2Vobj.edgeData(patch).v_xi{midx};
                    case {3,4}
                        v_tm = sgn_m*S2Vobj.edgeData(patch).v_eta{midx};
                end
                v_nm = S2Vobj.edgeData(patch).v_n{midx};
                sgn_s = (-1)^sidx;
                switch sidx
                    case {1,2}
                        v_ts = sgn_s*S2Vobj.edgeData(slave).v_xi{sidx};
                    case {3,4}
                        v_ts = sgn_s*S2Vobj.edgeData(slave).v_eta{sidx};
                end
                v_ns = S2Vobj.edgeData(slave).v_n{sidx};
                if orient
                    v_ts = flipud(v_ts);
                    v_ns = flipud(v_ns);
                end
            end
            nanIndices = or(any(isnan(v_tm),2), any(isnan(v_ts),2));
            v_tm(nanIndices,:) = [];
            v_ts(nanIndices,:) = [];
            v_nm(nanIndices,:) = [];
            v_ns(nanIndices,:) = [];

            angle = abs(acos(dot(v_nm,v_ns,2)));
            indices = dot((v_ts+v_tm)/2,(v_ns+v_nm)/2,2) < 0; % find edges forming a corner
            angle(indices) = pi - angle(indices);
            angle(~indices) = pi + angle(~indices); % the boundary at these corners are convex
            
            S2Vobj.angles(patch,midx) = max(angle);
            S2Vobj.angles(slave,sidx) = S2Vobj.angles(patch,midx);
        end
    end
end

function midx_regular = singularFromRegular(topologyMap,patch,midx)
% Check if singularity originates from a regular patch-grid
midx_regular = [];
for i = 1:numel(midx)
    if isempty(topologyMap{patch}(midx(i)).slave)
        neighbours = midx2neighbours(midx(i));
        slave1 = topologyMap{patch}(neighbours(1)).slave;
        sidx1 = topologyMap{patch}(neighbours(1)).sidx;
        slave2 = topologyMap{patch}(neighbours(2)).slave;
        sidx2 = topologyMap{patch}(neighbours(2)).sidx;
        counter1 = 0;
        while ~isempty(slave1) && slave1 ~= patch
            prev_slave1 = slave1;
            slave1 = topologyMap{slave1}(midxToRightmidx(sidx1)).slave;
            sidx1 = topologyMap{prev_slave1}(midxToRightmidx(sidx1)).sidx;
            counter1 = counter1 + 1;
        end
        counter2 = 0;
        while ~isempty(slave2) && slave2 ~= patch
            prev_slave2 = slave2;
            slave2 = topologyMap{slave2}(midxToLeftmidx(sidx2)).slave;
            sidx2 = topologyMap{prev_slave2}(midxToLeftmidx(sidx2)).sidx;
            counter2 = counter2 + 1;
        end
    
        if prev_slave1 == prev_slave2 && ~((~isempty(slave1) && slave1 == patch) || (~isempty(slave2) && slave2 == patch)) && (counter1 == 2 || counter2 == 2) % The final condition ensures that the new patch is connected to one of the slaves
            midx_regular = [midx_regular,midx(i)];
        end
    end
end

function midx = correct_midx_order(midx)
switch numel(midx)
    case 2
        if all(midx == [2,3])
            midx = [3,2];
        elseif all(midx == [1,4])
            midx = [4,1];
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
end

function i_corner = intermediateCorner(midx)

if all(ismember([1,3],midx))
    i_corner = 1;
elseif all(ismember([3,2],midx))
    i_corner = 2;
elseif all(ismember([2,4],midx))
    i_corner = 3;
else
    i_corner = 4;
end

function idx_oposite = get_oposite_corner_index(idx)
idx_oposite = mod(idx+2-1,4)+1;

function i = midx2leftCorner(midx)
map = [4,2,1,3]; % This map gives the i_corner for the corner (to the left) of the patch from a given side. 
i = map(midx);

function i = midx2rightCorner(midx)
map = [1,3,2,4]; % This map gives the i_corner for the corner (to the right) of the patch from a given side. 
i = map(midx);

function midx = midxToLeftmidx(sidx)
map = [4,3,1,2]; % This map gives the side index to the left of a given side
midx = map(sidx);

function midx = midxToRightmidx(sidx)
map = [3,4,2,1]; % This map gives the side index to the right of a given side
midx = map(sidx);

function midices = midx2neighbours(sidx)
map = [3,4; 4,3; 2,1; 1,2];
midices = map(sidx,:);

function indices = midx2corners(midx)
map = [1,4; 2,3; 1,2; 4,3];
indices = map(midx,:);

function midx = corners2midxes(i)
map = [1,3; 3,2; 2,4; 4,1];
midx = map(i,:);

function idx = corner2leftmidx(i)
map =  [1,3,2,4]; % This map gives the midx for the side of the patch sharing a given corner to the left of that corner. 
idx = map(i);

function idx = idx1To_idx2(idx1,idx2)
map = [1,3,2,4];
i_map1 = find(map == idx1);
map = circshift(map,1-i_map1);
idx = map([1,3,2,4] == idx2);

function orient = idx1To_idx2_orient(idx1,idx2,flip)
if nargin < 3
    flip = false;
end
if flip
    switch idx2
        case 1
            map =  [2,1,4,7]; 
        case 3
            map =  [4,7,1,2]; 
        otherwise
            error('not implemented')
    end
else
    switch idx2
        case 1
            map =  [0,3,5,6]; 
        case 3
            map =  [6,5,0,3]; 
        otherwise
            error('not implemented')
    end
end
orient = map(idx1);

function face = faceFromNormals(patch,S2Vobj,options,t,useAvgNormals)
if nargin < 5
    useAvgNormals = true;
end
nurbs = S2Vobj.nurbs_bdry{patch};
d = nurbs.d;
knots = nurbs.knots;
degree = nurbs.degree;
gxi = aveknt(knots{1}, degree(1)+1);
geta = aveknt(knots{2}, degree(2)+1);
[XI,ETA] = ndgrid(gxi,geta);
[~,dXdxi,dXdeta] = evaluateNURBSvec(nurbs,[XI(:),ETA(:)],1);
v_xi = dXdxi./norm2(dXdxi);
v_eta = dXdeta./norm2(dXdeta);
v_n = cross(v_xi, v_eta, 2);
v_n = v_n./norm2(v_n);
v_n = reshape(v_n,numel(gxi),numel(geta),d);

% Use corrected normal vectors at the boundaries
v_n(1,:,:) = S2Vobj.edgeData(patch).v_n{1};
v_n(end,:,:) = S2Vobj.edgeData(patch).v_n{2};
v_n(:,1,:) = S2Vobj.edgeData(patch).v_n{3};
v_n(:,end,:) = S2Vobj.edgeData(patch).v_n{4};

if useAvgNormals
    for midx2 = 1:4
        avg_v_n = S2Vobj.edgeData(patch).avg_v_n{midx2};
        v_n_temp = S2Vobj.edgeData(patch).v_n{midx2};
        indices = acos(abs(dot(avg_v_n,v_n_temp,2))) > options.avg_v_n_threshholdAngle;
        v_n_temp(indices,:) = v_n_temp(indices,:);
        switch midx2
            case 1
                v_n(1,:,:) = v_n_temp;
            case 2
                v_n(end,:,:) = v_n_temp;
            case 3
                v_n(:,1,:) = v_n_temp;
            case 4
                v_n(:,end,:) = v_n_temp;
        end
    end
end
coeffs = nurbs.coeffs;
if nargin < 4
    coeffs(1:3,:,:) = permute(v_n,[3,1,2]);
else
    coeffs(1:3,:,:) = coeffs(1:3,:,:) + t.*permute(v_n,[3,1,2]);
end
face = createNURBSobject(coeffs,knots);












