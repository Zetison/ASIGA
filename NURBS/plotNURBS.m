function [handles,maxC,minC] = plotNURBS(varargin)

%% Interpret input arguments
nurbsPatches = varargin{1};
if isempty(nurbsPatches)
    warning('No NURBS patches to be plotted')
    return
end
if isgraphics(nurbsPatches)
    ax = nurbsPatches;
    nurbsPatches = varargin{2};
    extraArg = 1;
    hold(ax,'on')
else
    ax = gca;
    extraArg = 0;
    hold on
end
% set default values
options = struct('resolution',[32,32,32], ...
                 'plotObject',true, ...
                 'plotElementEdges',true,...
                 'plotControlPolygon', false, ...
                 'plotParmDir', false, ...
                 'plotNormalVectors', false, ...
                 'plotColorFun', true, ...
                 'plotWeights', false, ...
                 'plotJacobian', false, ...
                 'coarseLinearSampling', true, ...
                 'color',jet(numel(nurbsPatches)),...
                 'alphaValue',1,...
                 'plotAt',true(3,2),...
                 'lineColor', 'black', ...
                 'quiverScale', NaN, ...
                 'quiverLineWidth', 2, ...
                 'qstep',NaN, ...
                 'quiverAutoScale','off', ...
                 'LineWidth',0.5, ...
                 'view', NaN, ...
                 'displayName', '', ...
                 'colorFun', NaN,...
                 'colorControlPolygon', 'red', ...
                 'markerEdgeColor', 'black', ...
                 'markerColor', 'black', ...
                 'UserData', [], ...
                 'app', [], ...
                 'colorParmDirs', {{getColor(12),getColor(13),getColor(7)}}, ...
                 'colorNormals', getColor(6));
if nargin > 1+extraArg
    if numel(varargin) > 2+extraArg
        newOptions = varargin(2+extraArg:end);
    else
        newOptions = varargin{2+extraArg};
    end
    options = updateOptions(options,newOptions);
end
if ~iscell(nurbsPatches)
    nurbsPatches = {nurbsPatches};
end
resolution = options.resolution;
plotObject = options.plotObject;
plotElementEdges = options.plotElementEdges;
plotControlPolygon = options.plotControlPolygon;
plotNormalVectors = options.plotNormalVectors;
plotParmDir = options.plotParmDir;
plotColorFun = options.plotColorFun;
quiverScale = options.quiverScale;
quiverAutoScale = options.quiverAutoScale;
if ~isa(quiverAutoScale,'char')
    if quiverAutoScale
        quiverAutoScale = 'on';
    else
        quiverAutoScale = 'off';
    end
end
if isempty(options.app)
    options.app.S2Vmouse_click = @(src,event) NaN;
end
if strcmp(quiverAutoScale,'off')
    quiverScaleAuto = 0.9;
else
    quiverScaleAuto = quiverScale;
    quiverScale = 1;
end
quiverLineWidth = options.quiverLineWidth;
markerEdgeColor = options.markerEdgeColor;
plotAt = options.plotAt;
colorFun = options.colorFun;
if ~isa(colorFun, 'function_handle')
    plotColorFun = false;
end
lineWidth = options.LineWidth;
lineColor = options.lineColor;
qstep = options.qstep;
colorControlPolygon = options.colorControlPolygon;
markerColor = options.markerColor;
color = options.color;
alphaValue = options.alphaValue;
colorParmDirs = options.colorParmDirs;
colorNormals = options.colorNormals;
if alphaValue ~= 1
    faceLighting = 'none';
else
    faceLighting = 'flat';
end
noPatches = numel(nurbsPatches);
if size(color,1) == 1 && noPatches > 1
    color = repmat(color,noPatches,1);
elseif size(color,1) < noPatches
    color = repmat(color,ceil(noPatches/size(color,1)),1);
end
plotJacobian = options.plotJacobian;
plotWeights = options.plotWeights;
if plotColorFun + plotJacobian + plotWeights > 1
    warning('Only one of the options plotColorFun, plotJacobian and plotWeights will be used')
end
% evaluateDerivatives = double(plotParmDir || plotNormalVectors || plotJacobian || (plotColorFun && nargin(colorFun) == 2));
evaluateDerivatives = 1;

visible = 'on'; % render only in the end, but it seems faster to render on the go, so this is here set to 'on'
maxC = NaN;
minC = NaN;
d_max = -Inf;
objectHandle = gobjects(1,noPatches);
elementEdgesHandle = gobjects(1,noPatches);
controlPolygonHandle = gobjects(1,noPatches);
parmDirHandle = gobjects(3,noPatches);
normalVectorsHandle = gobjects(1,noPatches);
for patch = 1:noPatches
    if isempty(options.displayName)
        displayName = sprintf('Patch %d ', patch);
    elseif isa(options.displayName, 'function_handle')
        displayName = options.displayName(patch);
    else
        displayName = [options.displayName ', patch ' num2str(patch)];
    end
    nurbs = nurbsPatches{patch};
    if isempty(nurbs)
        warning('Empty NURBS patch encountered')
        continue
    end
    if isfield(nurbs,'isPML') && any(nurbs.isPML)
        colorPatch = getColor(11);
    else
        colorPatch = color(patch,:);
    end
    d_p = nurbs.d_p;
    d = nurbs.d;
    if d > d_max
        d_max = d;
    end
    if plotColorFun && ((plotJacobian && d_p == 3) || plotWeights)
        warning('Plotting solution from function handle instead of Jacobian')
    end
    p_values = cell(1,d_p);
    uniqueKnots = cell(1,d_p);
    noElemsDir = zeros(1,d_p);
    noUniqueKnots = zeros(1,d_p);
    res = zeros(1,d_p);
    for i = 1:d_p
        % To reconstruct any interpolatory points, any p repeated knot should be
        % included in the abscissa values to plot
        uniqueKnots{i} = unique(nurbs.knots{i});
        noUniqueKnots(i) = numel(uniqueKnots{i});
        noElemsDir(i) = noUniqueKnots(i)-1;

        p = nurbs.degree(i);
        if p == 1 && options.coarseLinearSampling
            res(i) = 0;
        else
            noC0elements = numel(find(knotRepetitions(nurbs.knots{i},uniqueKnots{i}) >= nurbs.degree(i)))-1;
            res(i) = round(noC0elements*resolution(i)/noElemsDir(i));
        end
        p_values{i} = insertUniform(uniqueKnots{i}, res(i));
    end
    if isnan(qstep)
        qstep = max(round(max(res)/2),1);
    end
    if d_p == 3
        indicesMat = [1,2,3;
                      2,3,1;
                      3,1,2];
        indicesMat2 = [1,2,3;
                       3,1,2;
                       2,3,1];
    else
        indicesMat = [1,2;
                      2,1];
        indicesMat2 = [1,2;
                       2,1];
    end
    if d_p == 3
        nuk1 = length(p_values{1});
        nuk2 = length(p_values{2});
        nuk3 = length(p_values{3});
        maxnuk = max([nuk1,nuk2,nuk3]);
        vElementEdges = NaN(6*(maxnuk+1)*max(noUniqueKnots)*2,3);
        X = NaN(maxnuk,6*(maxnuk+1),3);
        normals = NaN(maxnuk,6*(maxnuk+1),3);
        Up = cell(3,3);
        for i_c = 1:d_p
            for j_c = 1:d_p
                Up{i_c,j_c} = NaN(maxnuk,6*(maxnuk+1));
            end
        end
        L_gamma = -Inf;
        C = NaN(maxnuk,6*(maxnuk+1));
        counter = 0;
        elCounter = 1;
        for ii = 1:d_p
            indices = indicesMat(ii,:);
            for jj = 1:size(plotAt,2)
                if plotAt(ii,jj)
                    nuk1 = length(p_values{indices(2)});
                    nuk2 = length(p_values{indices(3)});
                    XI = (jj-1)*ones(nuk1,nuk2);
                    [ETA,ZETA] = ndgrid(p_values{indices(2)},p_values{indices(3)});
                    XIETAZETA = [XI(:), ETA(:), ZETA(:)];
                    [X_temp,dvdxi,dvdeta,dvdzeta] = evaluateNURBS(nurbs, XIETAZETA(:,indicesMat2(ii,:)), evaluateDerivatives);
                    if evaluateDerivatives
                        dX_temp = zeros(nuk1*nuk2, d, d_p);
                        dX_temp(:,:,1) = dvdxi./norm2(dvdxi);
                        dX_temp(:,:,2) = dvdeta./norm2(dvdeta);
                        dX_temp(:,:,3) = dvdzeta./norm2(dvdzeta);
                        if plotNormalVectors || (plotColorFun && nargin(colorFun) == 2)
                            normal = cross(dX_temp(:,:,indices(2)),dX_temp(:,:,indices(3)),2);
                            normals_temp = (-1)^jj*normal./norm2(normal);
                        end
                        dX_temp = reshape(dX_temp,nuk1,nuk2, d, d_p);
                        if plotParmDir
                            for i_c = 1:d_p
                                for j_c = 1:d_p
                                    Up{i_c,j_c}(1:nuk1,counter+1:counter+nuk2) = dX_temp(:,:,i_c,j_c);
                                end
                            end
                        end
                        if plotJacobian
                            J_1 = dot(dX_temp(:,:,:,1),cross(dX_temp(:,:,:,2),dX_temp(:,:,:,3),3),3);
                            J_1(J_1 < eps) = eps;
                            if any(J_1(:) < -10*eps)
                                warning('Negative Jacobian encountered')
                            end
                            C(1:nuk1,counter+1:counter+nuk2) = J_1;
                        end
                    end
                    if plotWeights
                        C_temp = evaluateNURBS(nurbs, XIETAZETA(:,indicesMat2(ii,:)), 0, nurbs.coeffs(4,:).');
                        C_temp = reshape(C_temp,nuk1,nuk2);
                        C(1:nuk1,counter+1:counter+nuk2) = C_temp;
                    end
                    if plotColorFun
                        switch nargin(colorFun)
                            case 1
                                C_temp = colorFun(X_temp);
                            case 2
                                C_temp = colorFun(X_temp,normals_temp);
                            case 3
                                n = reshape(cross(dX_temp(:,:,:,1),dX_temp(:,:,:,2),3),nuk1*nuk2,3);
                                n = n./norm2(n);
                                n(isnan(n)) = 10;
                                dXdzeta = reshape(dX_temp(:,:,:,3),nuk1*nuk2,3);
                                C_temp = colorFun(X_temp,n,dXdzeta./norm2(dXdzeta));
                            case 4
                                C_temp = colorFun(XIETAZETA(:,indicesMat2(ii,:)),nurbs,NaN,NaN);
                        end
                        C(1:nuk1,counter+1:counter+nuk2) = reshape(C_temp,nuk1,nuk2);
                    end
                    L_gamma_temp = norm(X_temp(end,:)-X_temp(1,:));
                    if L_gamma_temp > L_gamma
                        L_gamma = L_gamma_temp;
                    end
                    X_temp = reshape(X_temp,nuk1,nuk2,d);

                    X(1:nuk1,counter+1:counter+nuk2,:) = X_temp;
                    if plotNormalVectors
                        normals(1:nuk1,counter+1:counter+nuk2,:) = reshape(normals_temp,nuk1,nuk2,d);
                    end
                    
                    % Create element edges
                    nuk3 = length(p_values{indices(1)});
                    uniqueKnots2 = nurbs.knots{indices(2)};
                    uniqueKnots3 = nurbs.knots{indices(3)};
                    noKnots2 = numel(uniqueKnots2);
                    noKnots3 = numel(uniqueKnots3);
                    [XI,ETA,ZETA] = ndgrid(p_values{indices(1)},uniqueKnots2,uniqueKnots3);
                    XIETAZETA = [XI(:), ETA(:), ZETA(:)];
                    X_temp2 = evaluateNURBS(nurbs, XIETAZETA(:,indicesMat2(ii,:)));
                    X_temp2 = reshape(X_temp2,nuk3,noKnots2,noKnots3,d);
                    X_temp2(end+1,:,:,:) = NaN;
                    noAddedPoints = (nuk3+1)*noKnots2*noKnots3;
                    vElementEdges(elCounter:elCounter+noAddedPoints-1,:) = reshape(X_temp2,[],3);
                    elCounter = elCounter+noAddedPoints;
                    
                    counter = counter + nuk2+1;
                end
            end
        end
        vElementEdges(elCounter:end,:) = [];

        if plotObject
            maxC = max(max(C));
            minC = min(min(C));
            if plotColorFun || plotJacobian || plotWeights
                objectHandle(patch) = surface(ax,X(:,:,1),X(:,:,2),X(:,:,3),C,'EdgeColor','none',...
                        'LineStyle','none','FaceAlpha',alphaValue, 'DisplayName',displayName,...
                        'Visible',visible,'ButtonDownFcn',@(src,event)options.app.S2Vmouse_click(src,event));
            else
                objectHandle(patch) = surface(ax,X(:,:,1),X(:,:,2),X(:,:,3), 'FaceColor', colorPatch,...
                        'EdgeColor','none','LineStyle','none','FaceAlpha',alphaValue, ...
                        'FaceLighting', faceLighting, 'DisplayName',displayName,...
                        'Visible',visible,'ButtonDownFcn',@(src,event)options.app.S2Vmouse_click(src,event));
            end
        end
        if plotElementEdges
            elementEdgesHandle(patch) = plotGridLines(vElementEdges,displayName);
        end
        if isnan(quiverScale)
            quiverScale = L_gamma/20;
        end
        if plotParmDir
            for i = 1:d_p
                if isfield(nurbs,'isPML') && nurbs.isPML(i)
                    colorParm = [1,0,0];
                else
                    colorParm = colorParmDirs{i};
                end
                parmDirHandle(i,patch) = quiver3(ax,X(1:qstep:end,1:qstep:end,1), ...
                        X(1:qstep:end,1:qstep:end,2), ...
                        X(1:qstep:end,1:qstep:end,3), ...
                        quiverScale*Up{1,i}(1:qstep:end,1:qstep:end), ...
                        quiverScale*Up{2,i}(1:qstep:end,1:qstep:end), ...
                        quiverScale*Up{3,i}(1:qstep:end,1:qstep:end),'color',colorParm,...
                        'LineWidth',quiverLineWidth,'AutoScale',quiverAutoScale,'AutoScaleFactor',quiverScaleAuto,...
                        'DisplayName',[displayName, ' - parm dir ' num2str(i)],'Visible',visible,'ButtonDownFcn',@(src,event)options.app.S2Vmouse_click(src,event));
            end
        end
        if plotNormalVectors
            normalVectorsHandle(patch) = quiver3(ax,X(1:qstep:end,1:qstep:end,1), ...
                    X(1:qstep:end,1:qstep:end,2), ...
                    X(1:qstep:end,1:qstep:end,3), ...
                    quiverScale*normals(1:qstep:end,1:qstep:end,1), ...
                    quiverScale*normals(1:qstep:end,1:qstep:end,2), ...
                    quiverScale*normals(1:qstep:end,1:qstep:end,3),'AutoScale',quiverAutoScale,'AutoScaleFactor',quiverScaleAuto,...
                    'DisplayName',[displayName, ' - normal vectors'],'color',colorNormals,'Visible',visible,'ButtonDownFcn',@(src,event)options.app.S2Vmouse_click(src,event));
        end
    elseif d_p == 2 && d == 3
        nuk1 = length(p_values{1});
        nuk2 = length(p_values{2});
        normals = zeros(length(p_values{1}), length(p_values{2}), d);
        C = zeros(length(p_values{1}), length(p_values{2}));
        [XI,ETA] = ndgrid(p_values{1},p_values{2});
        [X,dvdxi,dvdeta] = evaluateNURBS(nurbs, [XI(:) ETA(:)], evaluateDerivatives);
        L_gamma = norm(X(end,:)-X(1,:));
        if evaluateDerivatives
            dX = zeros(nuk1*nuk2, d, d_p);
            dX(:,:,1) = dvdxi./norm2(dvdxi);
            dX(:,:,2) = dvdeta./norm2(dvdeta);
            dX = reshape(dX,nuk1,nuk2, d, d_p);
            if plotNormalVectors || (plotColorFun && nargin(colorFun) == 2)
                normal = cross(dvdxi,dvdeta,2);
                normals = normal./norm2(normal);
            end
        end
        if plotWeights
            C = evaluateNURBS(nurbs, [XI(:) ETA(:)], 1, nurbs.coeffs(4,:).');
            C = reshape(C,nuk1,nuk2);
        end
        if plotColorFun
            switch nargin(colorFun)
                case 1
                    C = colorFun(X);
                case 2
                    C = colorFun(X,normals);
                case 4
                    C = colorFun([XI(:) ETA(:)],nurbs,NaN,NaN);
            end
            C = reshape(C,nuk1,nuk2);
        end
        X = reshape(X,nuk1,nuk2,d);
    
        if plotObject
            if plotColorFun || plotWeights
                objectHandle(patch) = surface(ax,X(:,:,1),X(:,:,2),X(:,:,3),C,'EdgeColor','none',...
                            'LineStyle','none','FaceAlpha',alphaValue, ...
                            'FaceLighting', faceLighting, 'DisplayName',displayName,...
                            'Visible',visible,'ButtonDownFcn',@(src,event)options.app.S2Vmouse_click(src,event));
                maxC = max(max(C));
                minC = min(min(C));
            else
                objectHandle(patch) = surface(ax,X(:,:,1),X(:,:,2),X(:,:,3), 'FaceColor', colorPatch,...
                            'EdgeColor','none','LineStyle','none','FaceAlpha',alphaValue, ...
                            'FaceLighting', faceLighting, 'DisplayName',displayName,...
                            'Visible',visible,'ButtonDownFcn',@(src,event)options.app.S2Vmouse_click(src,event));
            end
        end
        if plotElementEdges
            stepLen = res+1;
            vElementEdges = NaN((nuk1+1)*noUniqueKnots(2)+(nuk1+1)*noUniqueKnots(2),3);

            elCounter = 1;
            temp = NaN(nuk1+1, noUniqueKnots(2), 3);
            temp(1:end-1,:,:) = X(:,1:stepLen(2):end,:);
            noAddedPoints = noUniqueKnots(2)*(nuk1+1);
            vElementEdges(elCounter:elCounter+noAddedPoints-1,:) = reshape(temp,[],3);
            elCounter = elCounter+noAddedPoints;

            temp = NaN(nuk2+1, noUniqueKnots(1), 3);
            temp(1:end-1,:,:) = permute(X(1:stepLen(1):end,:,:),[2,1,3]);
            noAddedPoints = noUniqueKnots(1)*(nuk2+1);
            vElementEdges(elCounter:elCounter+noAddedPoints-1,:) = reshape(temp,[],3);
            elementEdgesHandle(patch) = plotGridLines(vElementEdges,displayName);
        end
        if isnan(quiverScale)
            quiverScale = L_gamma/20;
        end
        if plotParmDir
            for i = 1:d_p
                parmDirHandle(i,patch) = quiver3(ax,X(1:qstep:end,1:qstep:end,1), ...
                        X(1:qstep:end,1:qstep:end,2), ...
                        X(1:qstep:end,1:qstep:end,3), ...
                        quiverScale*dX(1:qstep:end,1:qstep:end,1,i), ...
                        quiverScale*dX(1:qstep:end,1:qstep:end,2,i), ...
                        quiverScale*dX(1:qstep:end,1:qstep:end,3,i),'LineWidth',quiverLineWidth,...
                        'color',colorParmDirs{i},'AutoScale',quiverAutoScale,'AutoScaleFactor',quiverScaleAuto,...
                        'DisplayName',[displayName, ' - parm dir ' num2str(i)],'Visible',visible,'ButtonDownFcn',@(src,event)options.app.S2Vmouse_click(src,event));
            end
        end
        if plotNormalVectors
            normals = reshape(normals,nuk1,nuk2,d);
            normalVectorsHandle(patch) = quiver3(ax,X(1:qstep:end,1:qstep:end,1), ...
                    X(1:qstep:end,1:qstep:end,2), ...
                    X(1:qstep:end,1:qstep:end,3), ...
                    quiverScale*normals(1:qstep:end,1:qstep:end,1), ...
                    quiverScale*normals(1:qstep:end,1:qstep:end,2), ...
                    quiverScale*normals(1:qstep:end,1:qstep:end,3),'AutoScale',quiverAutoScale,'AutoScaleFactor',quiverScaleAuto,...
                    'DisplayName',[displayName, ' - normal vectors'],'color',colorNormals,'Visible',visible,'ButtonDownFcn',@(src,event)options.app.S2Vmouse_click(src,event));
        end
    elseif d_p == 2 && d == 2
        if plotColorFun || plotParmDir
            nuk1 = length(p_values{1});
            nuk2 = length(p_values{2});
            C = zeros(length(p_values{1}), length(p_values{2}));
            [XI,ETA] = ndgrid(p_values{1},p_values{2});
            [X,dvdxi,dvdeta] = evaluateNURBS(nurbs, [XI(:) ETA(:)], evaluateDerivatives);
            L_gamma = norm(X(end,:)-X(1,:));
            if isnan(quiverScale)
                quiverScale = L_gamma/20;
            end
            if evaluateDerivatives
                dX = zeros(nuk1*nuk2, d, d_p);
                dX(:,:,1) = dvdxi./norm2(dvdxi);
                dX(:,:,2) = dvdeta./norm2(dvdeta);
                dX = reshape(dX,nuk1,nuk2, d, d_p);
            end
            if plotWeights
                C = evaluateNURBS(nurbs, [XI(:) ETA(:)], 0, nurbs.coeffs(4,:).');
                C = reshape(C,nuk1,nuk2);
            end
            if plotColorFun
                switch nargin(colorFun)
                    case 1
                        C = colorFun(X);
                    case 4
                        C = colorFun([XI(:) ETA(:)],nurbs,NaN,NaN);
                end
                C = reshape(C,nuk1,nuk2);
            end
            X = reshape(X,nuk1,nuk2,d);

            maxC = max(max(C));
            minC = min(min(C));
            if plotObject
                objectHandle(patch) = surface(ax,X(:,:,1),X(:,:,2),C-maxC,C,'EdgeColor','none',...
                            'LineStyle','none','FaceAlpha',alphaValue, 'FaceLighting', faceLighting, ...
                            'DisplayName',displayName,'Visible',visible,'ButtonDownFcn',@(src,event)options.app.S2Vmouse_click(src,event));
            end
            if plotParmDir
                for i = 1:d_p
                    parmDirHandle(i,patch) = quiver(X(1:qstep:end,1:qstep:end,1), ...
                            X(1:qstep:end,1:qstep:end,2), ...
                            quiverScale*dX(1:qstep:end,1:qstep:end,1,i), ...
                            quiverScale*dX(1:qstep:end,1:qstep:end,2,i),'LineWidth',quiverLineWidth,...
                            'color',colorParmDirs{i},'AutoScale',quiverAutoScale,'AutoScaleFactor',quiverScaleAuto,...
                            'DisplayName',[displayName, ' - parm dir ' num2str(i)],'Visible',visible,'ButtonDownFcn',@(src,event)options.app.S2Vmouse_click(src,event));
                end
            end
        else
            % Find the curve counter-clockwise around the patch
            XI = zeros(2*length(p_values{1})+2*length(p_values{2})-4, d);
            counter = 1;
            
            % along eta = 0
            xi = p_values{1};
            npts = size(xi,1);
            XI(counter:counter+npts-1,:) = [xi, zeros(npts,1)];
            counter = counter + npts;
    
            % along xi = 1
            eta = p_values{2}(2:end);
            npts = size(eta,1);
            XI(counter:counter+npts-1,:) = [ones(npts,1), eta];
            counter = counter + npts;
    
            % along eta = 1
            xi = p_values{1}(end-1:-1:1);
            npts = size(xi,1);
            XI(counter:counter+npts-1,:) = [xi, ones(npts,1)];
            counter = counter + npts;
            
            % along xi = 0
            eta = p_values{2}(end-1:-1:2); % Skip the end-point as the fill-function will automatically connect the curve
            npts = size(eta,1);
            XI(counter:counter+npts-1,:) = [zeros(npts,1), eta];
    
            % evaluate NURBS
            v = evaluateNURBS(nurbs, XI);
    
            if plotObject
                objectHandle(patch) = fill(ax,v(:,1),v(:,2), colorPatch,'EdgeColor','none','LineStyle','none', 'DisplayName',displayName,'Visible',visible,'ButtonDownFcn',@(src,event)options.app.S2Vmouse_click(src,event));
            end
        end

        if plotElementEdges
            nok1 = length(p_values{1}); % number of knots in the xi direction
            nok2 = length(p_values{2}); % number of knots in the eta direction
            nuk1 = length(uniqueKnots{1}); % number of unique knots in the xi direction
            nuk2 = length(uniqueKnots{2}); % number of unique knots in the eta direction
            XI = [[kron(ones(nuk2,1),p_values{1}), kron(uniqueKnots{2}.',ones(nok1,1))];
                  [kron(uniqueKnots{1}.',ones(nok2,1)), kron(ones(nuk1,1),p_values{2})]];
            npts = nuk2*(nok1+1)+nuk1*(nok2+1) - 1;
            v = NaN(npts, d);
            indices = false(npts,1);
            indices_NaN = [(nok1+1):(nok1+1):(nuk2*(nok1+1)), ...
                          ((nok2+1):(nok2+1):(nuk1*(nok2+1))-1)+(nuk2*(nok1+1))];
            indices(indices_NaN) = true;
            v(~indices,:) = evaluateNURBS(nurbs, XI);

            if plotElementEdges
                elementEdgesHandle(patch) = plotGridLines(v,displayName);
            end
        end
    elseif d_p == 1 && plotObject
        C = evaluateNURBS(nurbs, p_values{1});
        switch d
            case 1
                objectHandle(patch) = line(ax,C,zeros(size(C)), 'color', colorPatch, 'DisplayName',displayName,'Visible',visible,'ButtonDownFcn',@(src,event)options.app.S2Vmouse_click(src,event));  
            case 2
                objectHandle(patch) = line(ax,C(:,1), C(:,2), 'color', colorPatch, 'DisplayName',displayName,'Visible',visible,'ButtonDownFcn',@(src,event)options.app.S2Vmouse_click(src,event));  
            case 3
                objectHandle(patch) = line(ax,C(:,1), C(:,2), C(:,3), 'color', colorPatch, 'DisplayName',displayName,'Visible',visible,'ButtonDownFcn',@(src,event)options.app.S2Vmouse_click(src,event)); 
        end  
        C = evaluateNURBS(nurbs, uniqueKnots{1}.');
        if plotElementEdges
            elementEdgesHandle(patch) = line(ax,C(:,1), C(:,2), 'Marker', 'x', 'LineStyle', 'none', 'DisplayName',[displayName, ' - elements'],'Visible',visible,'ButtonDownFcn',@(src,event)options.app.S2Vmouse_click(src,event));
        end
    end
    if plotControlPolygon
        coeffs = subasgnArr(nurbs.coeffs,[],size(nurbs.coeffs,1)); % Remove the weights
        dimensions = size(coeffs);
        number = nurbs.number;

        if d_p == 1
            v = coeffs;
        else
            if d_p == 2
                v = zeros(d,d_p*prod(number)+sum(number));
            elseif d_p == 3
                v = zeros(d,3*prod(number)+prod(number([1,2]))+prod(number([1,3]))+prod(number([2,3])));
            end
            counter = 1;
            for i = 1:d_p
                indices = 1:d_p+1;
                indices([i+1,2]) = [2,i+1];
                temp = permute(coeffs,indices);
                prd = prod(dimensions(indices(3:end))); % product of remaining dimensions
                temp = cat(2,reshape(temp,d,number(i),prd),NaN(d,1,prd));  % adding NaN to make the data discontinous
                temp = reshape(temp,d,(number(i)+1)*prd);
                npts = size(temp,2) + 1;
                v(:,counter:counter+npts-1) = [temp, NaN(d,1)];
                counter = counter + npts;
            end
        end
        switch d
            case 1 % Need to use plot instead of line here in order to get LineStyle 'o-'
                controlPolygonHandle(patch) = line(ax,v,zeros(size(v)),'color',colorControlPolygon,'Marker','o','MarkerFaceColor', markerColor, ...
                            'MarkerEdgeColor', markerEdgeColor,'DisplayName',[displayName, ' - control polygon'],'Visible',visible,'ButtonDownFcn',@(src,event)options.app.S2Vmouse_click(src,event));
            case 2
                controlPolygonHandle(patch) = line(ax,v(1,:),v(2,:),'color',colorControlPolygon,'Marker','o','MarkerFaceColor', markerColor, ...
                                    'MarkerEdgeColor', markerEdgeColor,'DisplayName',[displayName, ' - control polygon'],'Visible',visible,'ButtonDownFcn',@(src,event)options.app.S2Vmouse_click(src,event));
            case 3
                controlPolygonHandle(patch) = line(ax,v(1,:),v(2,:),v(3,:),'color',colorControlPolygon,'Marker','o','MarkerFaceColor', markerColor, ...
                                            'MarkerEdgeColor', markerEdgeColor,'DisplayName',[displayName, ' - control polygon'],'Visible',visible,'ButtonDownFcn',@(src,event)options.app.S2Vmouse_click(src,event));
        end
    end
    if plotObject && ~isempty(options.UserData)
        options.UserData.patch = options.UserData.indices(patch);
        objectHandle(patch).UserData = options.UserData;
    end
    if plotElementEdges && ~isempty(options.UserData)
        elementEdgesHandle(patch).UserData = options.UserData;
    end
    if plotControlPolygon && ~isempty(options.UserData)
        controlPolygonHandle(patch).UserData = options.UserData;
    end
    if plotParmDir && ~isempty(options.UserData)
        for i = 1:size(parmDirHandle,1)
            if ~isa(parmDirHandle(i,patch),'matlab.graphics.GraphicsPlaceholder') 
                parmDirHandle(i,patch).UserData = options.UserData;
            end
        end
    end
    if plotNormalVectors && ~isempty(options.UserData)
        normalVectorsHandle(patch).UserData = options.UserData;
    end
end
if strcmp(visible,'off')
    for patch = 1:noPatches
        if ~isa(objectHandle(1),'matlab.graphics.GraphicsPlaceholder') 
            set(objectHandle,'Visible','on')
        end
        if ~isa(elementEdgesHandle(1),'matlab.graphics.GraphicsPlaceholder') 
            set(elementEdgesHandle,'Visible','on')
        end
        if ~isa(controlPolygonHandle(1),'matlab.graphics.GraphicsPlaceholder') 
            set(controlPolygonHandle,'Visible','on')
        end
        if ~isa(normalVectorsHandle(1),'matlab.graphics.GraphicsPlaceholder') 
            set(normalVectorsHandle,'Visible','on')
        end
        if ~isa(parmDirHandle(1),'matlab.graphics.GraphicsPlaceholder') 
            set(parmDirHandle,'Visible','on')
        end
    end
end
if extraArg == 1
    if d_max < 3
        view(ax,[0,90])
        axis(ax,'on')
    else
        axis(ax,'off')
        ax.Clipping = 'off';    % turn clipping off
    %     ax.SortMethod = 'ChildOrder';
        camproj(ax,'perspective')
        if ~isnan(options.view)
            view(ax,options.view)
        end
    end
    handles.ax = ax;
    if plotColorFun
        colormap(ax, 'parula');
    end
else
    if d_max < 3
        view(0,90)
        axis on
    else
        axis off
        ax = gca;               % get the current axis
        ax.Clipping = 'off';    % turn clipping off
    %     ax.SortMethod = 'ChildOrder';
        camproj('perspective')
        if ~isnan(options.view)
            view(options.view)
        end
    end
    handles.ax = gcf;
    if plotColorFun
        colorbar
    end
end
handles.objectHandle = objectHandle;
handles.elementEdgesHandle = elementEdgesHandle;
handles.controlPolygonHandle = controlPolygonHandle;
handles.parmDirHandle = parmDirHandle;
handles.normalVectorsHandle = normalVectorsHandle;
% axis equal

function h = plotGridLines(v,displayName)
if size(v,2) > 2
    h = line(ax,v(:,1),v(:,2),v(:,3),'color',lineColor,'LineWidth',lineWidth,'DisplayName',[displayName, ' - element edges'],'Visible',visible,'ButtonDownFcn',@(src,event)options.app.S2Vmouse_click(src,event));
else
    h = line(ax,v(:,1),v(:,2),'color',lineColor,'LineWidth',lineWidth,'DisplayName',[displayName, ' - element edges'],'Visible',visible,'ButtonDownFcn',@(src,event)options.app.S2Vmouse_click(src,event));
end
end
end
