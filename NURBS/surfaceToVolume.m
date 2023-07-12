function nurbs_vol = surfaceToVolume(varargin)
% This routine assumes nurbs_surf to be a closed NURBS-surface with all
% patches having normal vectors pointing outwards.

% set default values
options = struct('t', 1, ...                                % thickness of patch created from a surface patch having all angles larger than "sharpAngle" w.r.t. neighbouring patches
                 'sharpAngle', 120*pi/180,...               % Threshold for a "sharp" angle
                 'Eps', 1e-10, ...                          % Threshold for assuming two physical points in the l2-norm to be identical
                 'explodeNURBSsurface', true, ...           % Explode patches having C0-elements
                 'S2V_algorithm',{{'A1_1','A13_1'}}, ...    % Specify the desired algorithm to be used. Here 'A1_1' mean ("A"lgorithm for when only side patch 1 is known and the first of the available algorithms)
                 'avg_v_n_threshholdAngle', 11.25*pi/180);  % Threshold angle deviation between the normal vectors v_n avg_v_n

nurbs_surf = varargin{1};
if nargin > 1
    if numel(varargin) > 2
        newOptions = varargin(2:end);
    else
        newOptions = varargin{2};
    end
    options = updateOptions(options,newOptions);
end
nurbs_surf = cleanNURBS(nurbs_surf,[],1e-6);

if options.explodeNURBSsurface
    nurbs_surf = explodeNURBS(nurbs_surf);
end
if isfield(options,'connection')
    connection = options.connection;
else
    geometry = getTopology(nurbs_surf);
    connection = geometry.topology.connection;
end
if isfield(options,'ax')
    ax = options.ax;
else
    ax = gca;
end

d_p = nurbs_surf{1}.d_p;

noBdryPatches = numel(nurbs_surf);

topologyMap = createTopologyMap(connection,noBdryPatches,d_p);

[cornerData, angles] = computeCornerData(nurbs_surf,topologyMap,options);

nurbs_bdry = nurbs_surf;


nurbs_vol = cell(1,4*numel(nurbs_surf));
counter = 1;
% colors = colormap('hsv');
% colors = colormap('jet');
% getColorMap
% step = max(4,size(colors,1)/(2*noBdryPatches));
% colors = colors(1:step:end,:);
% colormap(flipud(colors))
% no_colors = size(colors,1);
no_colors = max(8,noBdryPatches/10);
colors = turbo(no_colors);
minC = Inf;
maxC = -Inf;


while counter == 1 || (any(angles(:) < pi) && any(~isnan(angles(1:noBdryPatches,:)),'all'))
%     if counter == 12
%         keyboard
%     end
    [nurbs_vol(counter), nurbs_covered,nurbs_newBdry] = addPatch(nurbs_bdry,topologyMap,cornerData,angles,options);
    
    if ~isempty(nurbs_newBdry)
        % Remove faces with zero measure
        zeroMeasure = NURBShasZeroMeasure(nurbs_newBdry);
        nurbs_newBdry(zeroMeasure) = [];
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Plotting
    if true
        plotNURBSvec(ax,nurbs_vol(counter),'plotControlPolygon',0,'plotParmDir',0,'displayName',['Patch ' num2str(counter)], ...
            'color',colors(mod(counter-1,no_colors)+1,:)); % colors(mod(round(counter*no_colors/(2*noBdryPatches))-1,no_colors)+1,:)
    else
        [~,maxC_patch,minC_patch] = plotNURBSvec(ax,nurbs_vol(counter),'plotControlPolygon',0,'plotNormalVectors',0,...
                        'displayName',['Patch ' num2str(counter)],'colorFun',@(xi,nurbs,b,c) meanRatioJacobian(nurbs,xi));
        if minC_patch < minC
            minC = minC_patch;
        end
        if maxC_patch > maxC
            maxC = maxC_patch;
        end
        clim(ax,[minC,maxC])
    end
    drawnow
    pause(0.001)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(nurbs_newBdry)
        noNewBdryPatches = numel(nurbs_newBdry);
        newBdryPatches = (numel(nurbs_bdry)+1):(numel(nurbs_bdry)+noNewBdryPatches);
    
        % Add new patches to the global set of surface patches, nurbs_bdry
        nurbs_bdry = [nurbs_bdry, nurbs_newBdry];
    
        % Find the surrounding patches to the patches that has been covered
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
    
    
        % Update topologyMap
        local2global = [surroundingPatches,newBdryPatches];
        geometry = getTopology(nurbs_bdry(local2global));
        connectionSub = geometry.topology.connection;
    
        for patch_i = newBdryPatches
            topologyMap{patch_i} = struct('slave', cell(1, 4), ...
                                          'sidx', cell(1, 4), ...
                                          'orient', cell(1, 4));
        end
    
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
    
        try
            % Update angles
            [cornerData,angles] = computeCornerData(nurbs_bdry,topologyMap,options,newBdryPatches,cornerData,angles);
        catch ME
            warning(ME.message)
            nurbs_vol(counter) = [];
            counter = counter - 1;
        end
    end
    angles(nurbs_covered,:) = NaN;
%     if counter == 1
%         for patch = 1:numel(cornerData)
%             if ~all(isnan(angles(patch,:)))
%                 quiver3(cornerData(patch).X(:,1), ...
%                         cornerData(patch).X(:,2), ...
%                         cornerData(patch).X(:,3), ...
%                         cornerData(patch).avg_v_n(:,1), ...
%                         cornerData(patch).avg_v_n(:,2), ...
%                         cornerData(patch).avg_v_n(:,3),'LineWidth',1,'AutoScale','off','DisplayName',['Patch ' num2str(patch)])
%             end
%         end
%         keyboard
%     end

    counter = counter + 1;

end
nurbs_vol(counter:end) = [];


function [nurbs, nurbs_covered,nurbs_newBdry] = addPatch(nurbs_bdry,topologyMap,cornerData,angles,options,patch,midx,maxNoSharpAngles)
%     [patch,maxNoSharpAngles,minSum] = findNextPatch(original_angles,sharpAngle); % Check for sharp angles in at the initial nurbs_surf first
%     midx = find(original_angles(patch,:) < sharpAngle);
%     if isinf(minSum)
%         [patch,maxNoSharpAngles] = findNextPatch(angles,sharpAngle);
%         midx = find(angles(patch,:) < sharpAngle);
%     end
d = 3; % dimension
t = options.t;
% original_angles = angles;
sharpAngle = options.sharpAngle;
avg_v_n_threshholdAngle = options.avg_v_n_threshholdAngle;
acuteAngleAdjustment = 1.2;
if nargin < 6
    [patch,maxNoSharpAngles] = findNextPatch(angles,sharpAngle);
    midx = find(angles(patch,:) < sharpAngle);
end
faces = cell(1,6);
switch maxNoSharpAngles
    case 0
        faces(1) = nurbs_bdry(patch);
        switch options.S2V_algorithm{maxNoSharpAngles+1}
            case 'A1_1'
                faces(1) = nurbs_bdry(patch);
                faces(2) = faceFromNormals(nurbs_bdry{patch},patch,cornerData,t,options);
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
        corners_opposite = midx2corners(midx-(-1)^midx);
        v_n_1 = cornerData(patch).v_n(corners_opposite,:);
        avg_v_n_i = cornerData(patch).avg_v_n(corners_opposite,:);
        angle_v_n_diff_i = mean(acos(abs(dot(v_n_1,avg_v_n_i,2))));

        corners_opposite = midx2corners(sidx-(-1)^sidx);
        v_n_1 = cornerData(slave).v_n(corners_opposite,:);
        avg_v_n_i = cornerData(slave).avg_v_n(corners_opposite,:);
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
        
        midx_opposite = midx-(-1)^midx;
        switch options.S2V_algorithm{maxNoSharpAngles+1}
            case 'A13_1'
                try
                    oppositeCornerIndices = midx2corners(midx_opposite);
                    if midx_opposite == 1 || midx_opposite == 4
                        oppositeCornerIndices = fliplr(oppositeCornerIndices);
                    end
                    cornerNormal = cornerData(patch).avg_v_n(oppositeCornerIndices,:).';
                    v_n = cornerData(patch).v_n(oppositeCornerIndices,:).';
                    if 1
                        indices = acos(abs(dot(cornerNormal,v_n,1))) > avg_v_n_threshholdAngle;
                        cornerNormal(:, indices) = v_n(:,indices);
                    else
                        if all(acos(abs(dot(cornerNormal,v_n,1))) > avg_v_n_threshholdAngle)
                            cornerNormal = v_n;
                        end
                    end
        
                    acuteConnection = (angles(patch,midx) + angles(patch,midx_opposite)/2) < pi;
                    if true
                        if acuteConnection
                            error('acute angles found')
                        end
                    else
                        if acuteConnection
                            sgn = (-1)^midx;
                            switch midx
                                case {1,2}
                                    v_t = cornerData(patch).v_xi(oppositeCornerIndices,:);
                                case {3,4}
                                    v_t = cornerData(patch).v_eta(oppositeCornerIndices,:);
                            end
                            v_n = cornerData(patch).v_n(oppositeCornerIndices,:);
                            theta = angles(patch,midx)*acuteAngleAdjustment;
        
                            cornerNormal = (-sgn*v_t*cos(theta) + v_n*sin(theta)).';
                        end
                    end
        
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
        
                    edges(1) = subNURBS(faces(3),'at',[0,0;0,1]);
                    edges(2) = subNURBS(faces(4),'at',[0,0;0,1]);
                    g = reshape(aveknt(knots{2}, degree(2)+1),1,[]);
                    coeffs3 = edges{1}.coeffs(:,1).*(1-g) + edges{2}.coeffs(:,1).*g;
                    edges(3) = createNURBSobject(coeffs3,knots{2});
                    coeffs4 = edges{1}.coeffs(:,end).*(1-g) + edges{2}.coeffs(:,end).*g;
                    edges(4) = createNURBSobject(coeffs4,knots{2});
                    faces(2) = GordonHall(edges);
        
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
                    
                    nurbs = GordonHall(faces);
                    
                    if checkOrientation(nurbs, 10)
                        error('Self intersection!')
                    end
                catch
                    try
                        len1 = norm(slave1_coeffs(1:3,end,1,1) - slave1_coeffs(1:3,1,1,1));
                        len2 = norm(slave1_coeffs(1:3,end,1,end) - slave1_coeffs(1:3,1,1,end));
                        t = mean([len1,len2]);
                        faces(2) = faceFromNormals(nurbs_bdry{patch},patch,cornerData,t,options);
                        faces{2} = orientNURBS(faces{2}, idx1To_idx2_orient(midx,1)); % Make "midx = 1"
                        nurbs = loftNURBS({faces(1),faces(2)},1,1);
                        nurbs{1}.coeffs(:,:,1,:) = slave1_coeffs(:,:,1,:);
                    catch
                        coeffs = master_coeffs + slave1_coeffs - slave1_coeffs(:,1,1,:);
                        nurbs = createNURBSobject(coeffs,knots);
                    end
                end
            case 'A13_2'
                coeffs = master_coeffs + slave1_coeffs - slave1_coeffs(:,1,1,:);
                nurbs = createNURBSobject(coeffs,knots);
        end
        nurbs_covered = [patch,slave];
        nurbs_newBdry = subNURBS(nurbs,'at',[0,1;0,1;1,1],'outwardPointingNormals',true);
    case 2
        midx_regular = singularFromRegular(topologyMap,patch,setdiff(1:4,midx),nurbs_bdry);
        if ~isempty(midx_regular)
            [nurbs, nurbs_covered,nurbs_newBdry] = addPatch(nurbs_bdry,topologyMap,cornerData,angles,options,patch,sort([midx,midx_regular]),maxNoSharpAngles+1);
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
                midx_opposite = midx-(-1).^midx;
                cornerIndices = intermediateCorner(midx);
                oppositeCornerIndices = intermediateCorner(midx_opposite);

                cornerNormal = cornerData(patch).avg_v_n(oppositeCornerIndices,:).';
                v_n = cornerData(patch).v_n(oppositeCornerIndices,:).';
                indices = acos(abs(dot(cornerNormal,v_n,1))) > avg_v_n_threshholdAngle;
                cornerNormal(:, indices) = v_n(:,indices);
                acuteConnection = (angles(patch,midx) + angles(patch,midx_opposite)/2) < pi;
                if true
                    if any(acuteConnection)
                        error('acute angles found')
                    end
                else
                    if all(acuteConnection) % Make opposite sides "parallel"
                        cornerNormal = cornerData(patch).v_t2(cornerIndices,:).';
                    else
                        for i = 1:2
                            sgn = (-1)^midx(i);
                            if acuteConnection(i)
                                switch midx(i)
                                    case {1,2}
                                        v_t = cornerData(patch).v_xi(oppositeCornerIndices,:);
                                    case {3,4}
                                        v_t = cornerData(patch).v_eta(oppositeCornerIndices,:);
                                end
                                v_n = cornerData(patch).v_n(oppositeCornerIndices,:);
%                                 phi = acos(dot(v_t,v_n));
                                theta = angles(patch,midx(i))*acuteAngleAdjustment;
%                                 cornerNormal = ([v_t; v_n]\[cos(theta); cos(theta-phi)]);
%                                 cornerNormal = cornerNormal./vecnorm(cornerNormal,1);
        
                                cornerNormal = (-sgn*v_t*cos(theta) + v_n*sin(theta)).';
                            end
                        end
                    end
                end

                len2 = norm(slave1_coeffs(1:3,end,1,end) - slave1_coeffs(1:3,1,1,end));
                g = reshape(aveknt(knots{1}, degree(1)+1),1,[]);
                coeffs2 = zeros(d+1,number(1));
                coeffs2(1:3,:) = master_coeffs(1:3,1,end,end) + len2*cornerNormal.*g;
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
                
                nurbs = GordonHall(faces);
                if checkOrientation(nurbs, 10)
                    error('Self intersection!')
                end
            catch
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
            [nurbs, nurbs_covered,nurbs_newBdry] = addPatch(nurbs_bdry,topologyMap,cornerData,angles,options,patch,sort([midx,midx_regular]),maxNoSharpAngles+1);
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

function [patch,maxNoSharpAngles,minSum] = findNextPatch(angles,sharpAngle)

%     midx = find(angles(patch,:) < sharpAngle);

sharpAngles = angles < sharpAngle;
% midx = NaN(size(sharpAngles));
% midx(sharpAngles) = find(sharpAngles);
indices_op = or(and(sharpAngles(:,1),sharpAngles(:,2)),and(sharpAngles(:,3),sharpAngles(:,4)));
indices_nan = and(any(isnan(angles),2),~all(isnan(angles),2));

noSharpAngles = sum(sharpAngles,2);
noSharpAngles(indices_op) = noSharpAngles(indices_op) + 0.1; % Prioritize cases where oposite faces are known
noSharpAngles(indices_nan) = noSharpAngles(indices_nan) + 0.2; % Prioritize cases containing singularities
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
maxNoSharpAngles = floor(maxNoSharpAngles);
[minSum,patch] = min(avgAngles);

function [cornerData,angles] = computeCornerData(nurbs,topologyMap,options,indices,cornerData,angles)
noBdryPatches = numel(nurbs);
d_p = nurbs{1}.d_p;
d = 3; % dimension
if nargin < 5
    indices = 1:noBdryPatches;
end
if nargin < 6
    cornerData = struct('X', cell(1, noBdryPatches), ...
                        'v_t1', cell(1, noBdryPatches), ...
                        'v_t2', cell(1, noBdryPatches), ...
                        'v_xi', cell(1, noBdryPatches), ...
                        'v_eta', cell(1, noBdryPatches), ...
                        'cornerAngle', cell(1, noBdryPatches), ...
                        'v_n', cell(1, noBdryPatches), ...
                        'avg_v_n', cell(1, noBdryPatches));
    angles = NaN(noBdryPatches,2*d_p);
else
    noNewFields = noBdryPatches-numel(cornerData);
    cornerData = [cornerData, struct('X', cell(1, noNewFields), ...
                                     'v_t1', cell(1, noNewFields), ...
                                     'v_t2', cell(1, noNewFields), ...
                                     'v_xi', cell(1, noNewFields), ...
                                     'v_eta', cell(1, noNewFields), ...
                                     'cornerAngle', cell(1, noNewFields), ...
                                     'v_n', cell(1, noNewFields), ...
                                     'avg_v_n', cell(1, noNewFields))];
    angles = [angles; NaN(noNewFields,2*d_p)];
end
Eps = options.Eps;


xi = [0,0;
      1,0;
      1,1;
      0,1]; % Corner indices are oriented counter-clockwise starting at xi,eta=0

%% Compute normal vector at the vertices
for patch = indices
    [X,dXdxi,dXdeta] = evaluateNURBSvec(nurbs{patch},xi,1);
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
        [~,dXdxi_m2,dXdeta_m2] = evaluateNURBSvec(nurbs{patch},xi_eps(singularCorners,:),1);
        v_xi(singularCorners,:) = dXdxi_m2./norm2(dXdxi_m2);
        v_eta(singularCorners,:) = dXdeta_m2./norm2(dXdeta_m2);
        if ~isempty(singularCorners_v_n)
            v_n(singularCorners_v_n,:) = cross(v_xi(singularCorners_v_n,:), v_eta(singularCorners_v_n,:), 2);
        end
    end
    cornerAngle = acos(dot(v_t1,v_t2,2));

    v_n = v_n./norm2(v_n);
    cornerData(patch).X = X;
    cornerData(patch).v_xi = v_xi;
    cornerData(patch).v_eta = v_eta;
    cornerData(patch).v_t1 = v_t1;
    cornerData(patch).v_t2 = v_t2;
    cornerData(patch).v_n = v_n;
    cornerData(patch).cornerAngle = cornerAngle;
    if any(isnan(cornerData(patch).v_n(:)))
        error('Could not compute normal vector at corner. Case not accounted for')
    end
end

%% Compute average normal vector in all vertices
noCorners = 2^d_p;
for patch = indices
    cornerData(patch).avg_v_n = zeros(4,d);
    for i_corner = 1:noCorners
        v_n = zeros(1,d);
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

            if ~any(isnan(cornerData(slave).v_n))
                cornerAngle = cornerData(slave).cornerAngle(slave_i_corner);
                v_n = v_n + cornerAngle*cornerData(slave).v_n(slave_i_corner,:)/(2*pi); % weight the normal vector with the angle of the respetive corner
            end
            slave_prev = slave;
            midx = corner2midx(slave_i_corner);
            slave_temp = topologyMap{slave}(midx).slave;
            if isempty(slave_temp)
                idx = mod(find(corner2midx(1:4) == midx)-2,4)+1;
                midx = corner2midx(idx);
                slave = topologyMap{slave}(midx).slave;
            else
                slave = slave_temp;
            end
            sidx = topologyMap{slave_prev}(midx).sidx;

            slave_i_corner = sidx2leftCorner(sidx);
            neighbours(counter,1) = slave;
            neighbours(counter,2) = slave_i_corner;
        end
        neighbours(counter+1,1) = patch;
        neighbours(counter+1,2) = i_corner;
        v_n = v_n./norm2(v_n); % Normalize normal vector
        
        % Update the average normal vector for all patches sharing this corner
        for i = 1:counter+1
            slave = neighbours(i,1);
            slave_i_corner = neighbours(i,2);
            X = cornerData(slave).X(slave_i_corner,:);
            for j = 1:4
                if norm(X-cornerData(slave).X(j,:)) < Eps
                    cornerData(slave).avg_v_n(j,:) = v_n;
                end
            end
        end
    end
end
%% Compute angle between all pairs of connected patches
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
            i_corners_m = midx2corners(midx);
            switch midx
                case {1,2}
                    v_tm = sgn_m*cornerData(patch).v_xi(i_corners_m,:);
                case {3,4}
                    v_tm = sgn_m*cornerData(patch).v_eta(i_corners_m,:);
            end
            v_nm = cornerData(patch).v_n(i_corners_m,:);
            i_corners_s = midx2corners(sidx);
            sgn_s = (-1)^sidx;
            switch sidx
                case {1,2}
                    v_ts = sgn_s*cornerData(slave).v_xi(i_corners_s,:);
                case {3,4}
                    v_ts = sgn_s*cornerData(slave).v_eta(i_corners_s,:);
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

            angle = abs(acos(dot(v_nm,v_ns,2)));
            indices = dot((v_ts+v_tm)/2,(v_ns+v_nm)/2,2) < 0; % find edges forming a corner
            angle(indices) = pi - angle(indices);
            angle(~indices) = pi + angle(~indices); % the boundary at these corners are convex
            
            angles(patch,midx) = max(angle);
            angles(slave,sidx) = angles(patch,midx);
        end
    end
end
if 0
    for patch = 1:noBdryPatches
        if ~all(isnan(angles(patch,:)))
            quiver3(cornerData(patch).X(:,1), ...
                    cornerData(patch).X(:,2), ...
                    cornerData(patch).X(:,3), ...
                    cornerData(patch).avg_v_n(:,1), ...
                    cornerData(patch).avg_v_n(:,2), ...
                    cornerData(patch).avg_v_n(:,3),'LineWidth',1,'AutoScale','off','DisplayName',['Patch ' num2str(patch)])
        end
    end
    keyboard
end

function midx_regular = singularFromRegular(topologyMap,patch,midx,nurbs_bdry)
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

function i = sidx2leftCorner(sidx)
map = [4,2,1,3]; % This map gives the i_corner for the corner (to the left) of the patch from a given side. 
i = map(sidx);

function i = sidx2rightCorner(sidx)
map = [1,3,2,4]; % This map gives the i_corner for the corner (to the right) of the patch from a given side. 
i = map(sidx);

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

function idx = corner2midx(i)
map =  [1,3,2,4]; % This map gives the midx for the side of the patch sharing a given corner to the left of that side. 
idx = map(i);

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

function face = faceFromNormals(nurbs,patch,cornerData,t,options)
d = nurbs.d;
Eps = options.Eps;
knots = nurbs.knots;
degree = nurbs.degree;
gxi = aveknt(knots{1}, degree(1)+1);
geta = aveknt(knots{2}, degree(2)+1);
[XI,ETA] = ndgrid(gxi,geta);
[X,dXdxi,dXdeta] = evaluateNURBSvec(nurbs,[XI(:),ETA(:)],1);
v_xi = dXdxi./norm2(dXdxi);
v_eta = dXdeta./norm2(dXdeta);
v_n = cross(v_xi, v_eta, 2);
v_n = v_n./norm2(v_n);
v_n = reshape(v_n,numel(gxi),numel(geta),d);
X = reshape(X,numel(gxi),numel(geta),d);

% Handle edges collapsed into a single point
if norm(reshape(X(1,1,:)-X(end,1,:),d,[])) < Eps
    v_n(:,1,:) = repmat(cornerData(patch).v_n(1,:),numel(gxi),1);
end
if norm(reshape(X(end,1,:)-X(end,end,:),d,[])) < Eps
    v_n(end,:,:) = repmat(cornerData(patch).v_n(2,:),numel(geta),1);
end
if norm(reshape(X(1,end,:)-X(end,end,:),d,[])) < Eps
    v_n(:,end,:) = repmat(cornerData(patch).v_n(3,:),numel(gxi),1);
end
if norm(reshape(X(1,1,:)-X(1,end,:),d,[])) < Eps
    v_n(1,:,:) = repmat(cornerData(patch).v_n(4,:),numel(geta),1);
end
coeffs = nurbs.coeffs;
coeffs(1:3,:,:) = coeffs(1:3,:,:) + t*permute(v_n,[3,1,2]);
face = createNURBSobject(coeffs,knots);