function [h,maxC,minC] = plotNURBSvec(varargin)

%% Interpret input arguments
nurbsPatches = varargin{1};
% set default values
options = struct('resolution',[32,32,32], ...
                 'plotObject',true, ...
                 'plotElementEdges',true,...
                 'plotControlPolygon', 0, ...
                 'plotNormalVectors', 0, ...
                 'plotJacobian', 0, ...
                 'plotParmDir', 0, ...
                 'coarseLinearSampling', true, ...
                 'color',jet(numel(nurbsPatches)),...
                 'alphaValue',1,...
                 'plotAt',true(3,2),...
                 'colorFun', NaN,...
                 'lineColor', 'black', ...
                 'colorControlPolygon', 'red', ...
                 'markerEdgeColor', 'black', ...
                 'markerColor', 'black', ...
                 'quiverScale', NaN, ...
                 'quiverLineWidth', 2, ...
                 'LineWidth',0.5, ...
                 'qstep',NaN, ...
                 'view', NaN, ...
                 'displayName', '', ...
                 'colorParmDirs', {{getColor(12),getColor(13),getColor(7)}}, ...
                 'samplingDistance', NaN, ...
                 'elementBasedSamples',false);
if nargin > 1
    if numel(varargin) > 2
        newOptions = varargin(2:end);
    else
        newOptions = varargin{2};
    end
    options = updateOptions(options,newOptions);
end
if ~iscell(nurbsPatches)
    nurbsPatches = {nurbsPatches};
end
resolution = options.resolution;
plotElementEdges = options.plotElementEdges;
plotControlPolygon = options.plotControlPolygon;
plotNormalVectors = options.plotNormalVectors;
plotParmDir = options.plotParmDir;
quiverScale = options.quiverScale;
quiverLineWidth = options.quiverLineWidth;
markerEdgeColor = options.markerEdgeColor;
plotAt = options.plotAt;
colorFun = options.colorFun;
lineWidth = options.LineWidth;
lineColor = options.lineColor;
qstep = options.qstep;
colorControlPolygon = options.colorControlPolygon;
markerColor = options.markerColor;
color = options.color;
alphaValue = options.alphaValue;
colorParmDirs = options.colorParmDirs;
if alphaValue ~= 1
    faceLighting = 'none';
else
    faceLighting = 'flat';
end
noPatches = numel(nurbsPatches);
if size(color,1) == 1 && noPatches > 1
    color = repmat(color,noPatches,1);
elseif size(color,1) < noPatches
    error('A color must exist for each patch')
end
plotJacobian = options.plotJacobian;
plotSolution = isa(colorFun, 'function_handle');
maxC = NaN;
minC = NaN;

hold on
d_max = -Inf;
for patch = 1:noPatches
    if isempty(options.displayName)
        displayName = sprintf('Patch %d ', patch);
    else
        displayName = [options.displayName ', patch ' num2str(patch)];
    end
    nurbs = nurbsPatches{patch};
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
    if plotSolution && plotJacobian && d_p == 3
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
            res(i) = round(resolution(i)/noElemsDir(i));
        end
        p_values{i} = insertUniform(uniqueKnots{i}, res(i));
    end
    if isnan(qstep)
        qstep = max(round(res(1)/2),1);
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
                    [X_temp,dvdxi,dvdeta,dvdzeta] = evaluateNURBSvec(nurbs, XIETAZETA(:,indicesMat2(ii,:)), 1);
                    dX_temp = zeros(nuk1*nuk2, d, d_p);
                    dX_temp(:,:,1) = dvdxi./norm2(dvdxi);
                    dX_temp(:,:,2) = dvdeta./norm2(dvdeta);
                    dX_temp(:,:,3) = dvdzeta./norm2(dvdzeta);
                    if plotNormalVectors || (plotSolution && nargin(colorFun) == 2)
                        normal = cross(dX_temp(:,:,indices(2)),dX_temp(:,:,indices(3)),2);
                        normals_temp = (-1)^jj*normal./norm2(normal);
                    end
                    dX_temp = reshape(dX_temp,nuk1,nuk2, d, d_p);
                    if plotSolution
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
                        C(1:nuk1,counter+1:counter+nuk2) = log10(J_1);
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
                    
                    stepLen = res(indices(2:3))+1;
                    noLines = noUniqueKnots(indices(2:3));
                    
                    % Create element edges
                    temp = NaN(nuk1+1, noLines(2), 3);
                    temp(1:end-1,:,:) = X_temp(:,1:stepLen(2):end,:);
                    noAddedPoints = noLines(2)*(nuk1+1);
                    vElementEdges(elCounter:elCounter+noAddedPoints-1,:) = reshape(temp,[],3);
                    elCounter = elCounter+noAddedPoints;
                    
                    temp = NaN(nuk2+1, noLines(1), 3);
                    temp(1:end-1,:,:) = permute(X_temp(1:stepLen(1):end,:,:),[2,1,3]);
                    noAddedPoints = noLines(1)*(nuk2+1);
                    vElementEdges(elCounter:elCounter+noAddedPoints-1,:) = reshape(temp,[],3);
                    elCounter = elCounter+noAddedPoints;
                    
                    counter = counter + nuk2+1;
                end
            end
        end
        vElementEdges(elCounter:end,:) = [];

        if options.plotObject
            maxC = max(max(C));
            minC = min(min(C));
            if plotSolution || plotJacobian
                surf(X(:,:,1),X(:,:,2),X(:,:,3),C,'EdgeColor','none','LineStyle','none','FaceAlpha',alphaValue, 'DisplayName',displayName)
                colorbar
            else
                surf(X(:,:,1),X(:,:,2),X(:,:,3), 'FaceColor', colorPatch,'EdgeColor','none','LineStyle','none','FaceAlpha',alphaValue, ...
                        'FaceLighting', faceLighting, 'DisplayName',displayName)
            end
        end
        if plotElementEdges
            plotGridLines(vElementEdges,displayName)
        end
        if isnan(quiverScale)
            quiverScale = L_gamma/20;
        end
        if plotNormalVectors
            quiver3(X(1:qstep:end,1:qstep:end,1), ...
                    X(1:qstep:end,1:qstep:end,2), ...
                    X(1:qstep:end,1:qstep:end,3), ...
                    quiverScale*normals(1:qstep:end,1:qstep:end,1), ...
                    quiverScale*normals(1:qstep:end,1:qstep:end,2), ...
                    quiverScale*normals(1:qstep:end,1:qstep:end,3),'AutoScale','off','DisplayName',[displayName, ' - normal vectors'])
        end
        if plotParmDir
            for i = 1:d_p
                if isfield(nurbs,'isPML') && nurbs.isPML(i)
                    colorParm = [1,0,0];
                else
                    colorParm = colorParmDirs{i};
                end
                quiver3(X(1:qstep:end,1:qstep:end,1), ...
                        X(1:qstep:end,1:qstep:end,2), ...
                        X(1:qstep:end,1:qstep:end,3), ...
                        quiverScale*Up{1,i}(1:qstep:end,1:qstep:end), ...
                        quiverScale*Up{2,i}(1:qstep:end,1:qstep:end), ...
                        quiverScale*Up{3,i}(1:qstep:end,1:qstep:end),'color',colorParm,'LineWidth',quiverLineWidth,'AutoScale','off','DisplayName',[displayName, ' - parm dir ' num2str(i)])
            end
        end
    elseif d_p == 2 && d == 3
        nuk1 = length(p_values{1});
        nuk2 = length(p_values{2});
        normals = zeros(length(p_values{1}), length(p_values{2}), d);
        C = zeros(length(p_values{1}), length(p_values{2}));
        [XI,ETA] = ndgrid(p_values{1},p_values{2});
        [X,dvdxi,dvdeta] = evaluateNURBSvec(nurbs, [XI(:) ETA(:)], 1);
        L_gamma = norm(X(end,:)-X(1,:));
        dX = zeros(nuk1*nuk2, d, d_p);
        dX(:,:,1) = dvdxi./norm2(dvdxi);
        dX(:,:,2) = dvdeta./norm2(dvdeta);
        dX = reshape(dX,nuk1,nuk2, d, d_p);
        if plotNormalVectors || (plotSolution && nargin(colorFun) == 2)
            normal = cross(dvdxi,dvdeta,2);
            normals = normal./norm2(normal);
        end
        if plotSolution
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
    
        if options.plotObject
            if plotSolution
                surf(X(:,:,1),X(:,:,2),X(:,:,3),C,'EdgeColor','none','LineStyle','none','FaceAlpha',alphaValue, ...
                            'FaceLighting', faceLighting, 'DisplayName',displayName)
                maxC = max(max(C));
                minC = min(min(C));
                colorbar
            else
                surf(X(:,:,1),X(:,:,2),X(:,:,3), 'FaceColor', colorPatch,'EdgeColor','none','LineStyle','none','FaceAlpha',alphaValue, ...
                            'FaceLighting', faceLighting, 'DisplayName',displayName)
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
            plotGridLines(vElementEdges,displayName)
        end
        if isnan(quiverScale)
            quiverScale = L_gamma/20;
        end
        if plotNormalVectors
            normals = reshape(normals,nuk1,nuk2,d);
            quiver3(X(1:qstep:end,1:qstep:end,1), ...
                    X(1:qstep:end,1:qstep:end,2), ...
                    X(1:qstep:end,1:qstep:end,3), ...
                    quiverScale*normals(1:qstep:end,1:qstep:end,1), ...
                    quiverScale*normals(1:qstep:end,1:qstep:end,2), ...
                    quiverScale*normals(1:qstep:end,1:qstep:end,3),'AutoScale','off','DisplayName',[displayName, ' - normal vectors'])
        end
        if plotParmDir
            for i = 1:d_p
                quiver3(X(1:qstep:end,1:qstep:end,1), ...
                        X(1:qstep:end,1:qstep:end,2), ...
                        X(1:qstep:end,1:qstep:end,3), ...
                        quiverScale*dX(1:qstep:end,1:qstep:end,1,i), ...
                        quiverScale*dX(1:qstep:end,1:qstep:end,2,i), ...
                        quiverScale*dX(1:qstep:end,1:qstep:end,3,i),'LineWidth',quiverLineWidth,'color',colorParmDirs{i},'AutoScale','off','DisplayName',[displayName, ' - parm dir ' num2str(i)])
            end
        end
    elseif d_p == 2 && d == 2
        if options.plotObject
            if plotSolution || plotParmDir
                nuk1 = length(p_values{1});
                nuk2 = length(p_values{2});
                C = zeros(length(p_values{1}), length(p_values{2}));
                [XI,ETA] = ndgrid(p_values{1},p_values{2});
                [X,dvdxi,dvdeta] = evaluateNURBSvec(nurbs, [XI(:) ETA(:)], 1);
                L_gamma = norm(X(end,:)-X(1,:));
                if isnan(quiverScale)
                    quiverScale = L_gamma/20;
                end
                dX = zeros(nuk1*nuk2, d, d_p);
                dX(:,:,1) = dvdxi./norm2(dvdxi);
                dX(:,:,2) = dvdeta./norm2(dvdeta);
                dX = reshape(dX,nuk1,nuk2, d, d_p);
                if plotSolution
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
                surf(X(:,:,1),X(:,:,2),C-maxC,C,'EdgeColor','none','LineStyle','none','FaceAlpha',alphaValue, 'FaceLighting', faceLighting, 'DisplayName',displayName)
                colorbar
                if plotParmDir
                    for i = 1:d_p
                        quiver(X(1:qstep:end,1:qstep:end,1), ...
                                X(1:qstep:end,1:qstep:end,2), ...
                                quiverScale*dX(1:qstep:end,1:qstep:end,1,i), ...
                                quiverScale*dX(1:qstep:end,1:qstep:end,2,i),'LineWidth',quiverLineWidth,'color',colorParmDirs{i},'AutoScale','off','DisplayName',[displayName, ' - parm dir ' num2str(i)])
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
                v = evaluateNURBSvec(nurbs, XI);
        
                % reverse back order of arrays
                fill(v(:,1),v(:,2), colorPatch,'EdgeColor','none','LineStyle','none', 'DisplayName',displayName)
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
            v(~indices,:) = evaluateNURBSvec(nurbs, XI);

            if plotElementEdges
                plotGridLines(v,displayName)
            end
        end
    elseif d_p == 1 && options.plotObject
        C = evaluateNURBSvec(nurbs, p_values{1});
        switch d
            case 1
                plot(C,zeros(size(C)), 'color', colorPatch, 'DisplayName',displayName);  
            case 2
                plot(C(:,1), C(:,2), 'color', colorPatch, 'DisplayName',displayName);  
            case 3
                plot3(C(:,1), C(:,2), C(:,3), 'color', colorPatch, 'DisplayName',displayName); 
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
            case 1
                plot(v,zeros(size(v)),'o-','color',colorControlPolygon,'MarkerFaceColor', markerColor, ...
                            'MarkerEdgeColor', markerEdgeColor,'DisplayName',[displayName, ' - control polygon'])
            case 2
                plot(v(1,:),v(2,:),'o-','color',colorControlPolygon,'MarkerFaceColor', markerColor, ...
                                    'MarkerEdgeColor', markerEdgeColor,'DisplayName',[displayName, ' - control polygon'])
            case 3
                plot3(v(1,:),v(2,:),v(3,:),'o-','color',colorControlPolygon,'MarkerFaceColor', markerColor, ...
                                            'MarkerEdgeColor', markerEdgeColor,'DisplayName',[displayName, ' - control polygon'])
        end
    end
end
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
% axis equal
drawnow
h = gcf;

function plotGridLines(v,displayName)
if size(v,2) > 2
    plot3(v(:,1),v(:,2),v(:,3),'color',lineColor,'LineWidth',lineWidth,'DisplayName',[displayName, ' - element edges'])
else
    plot(v(:,1),v(:,2),'color',lineColor,'LineWidth',lineWidth,'DisplayName',[displayName, ' - element edges'])
end
end
end
