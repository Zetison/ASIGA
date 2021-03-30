function [h,maxC,minC] = plotNURBS(varargin)

%% Interpret input arguments
nurbsPatches = varargin{1};
% set default values
options = struct('resolution',[10,10,10], ...
                 'plotObject',true, ...
                 'plotElementEdges',true,...
                 'plotControlPolygon', 0, ...
                 'plotNormalVectors', 0, ...
                 'plotJacobian', 0, ...
                 'plotParmDir', 0, ...
                 'quiverScale', 0.2, ...
                 'coarseLinearSampling', true, ...
                 'color',jet(numel(nurbsPatches)),...
                 'alphaValue',1,...
                 'plotAt',true(3,2),...
                 'colorFun', NaN,...
                 'lineColor', 'black', ...
                 'colorControlPolygon', 'red', ...
                 'markerEdgeColor', 'black', ...
                 'markerColor', 'black', ...
                 'LineWidth',0.5, ...
                 'view', getView(), ...
                 'displayName', '', ...
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
markerEdgeColor = options.markerEdgeColor;
plotAt = options.plotAt;
colorFun = options.colorFun;
lineWidth = options.LineWidth;
lineColor = options.lineColor;
colorControlPolygon = options.colorControlPolygon;
markerColor = options.markerColor;
color = options.color;
alphaValue = options.alphaValue;
colorParmDirs = {getColor(12),getColor(13),getColor(7)};
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
    if isfield(nurbs,'isPML') && nurbs.isPML
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
    coeffs = nurbs.coeffs;
    coeffs = subasgnArr(coeffs,[],size(coeffs,1));
    dimensions = size(coeffs);
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
    if options.plotObject
        if d_p == 3
            nuk1 = length(p_values{1});
            nuk2 = length(p_values{2});
            nuk3 = length(p_values{3});
            maxnuk = max([nuk1,nuk2,nuk3]);
            vElementEdges = NaN(6*(maxnuk+1)*max(noUniqueKnots)*2,3);
            X = NaN(maxnuk,6*(maxnuk+1),3);
            Up = cell(3,3);
            for i_c = 1:d_p
                for j_c = 1:d_p
                    Up{i_c,j_c} = NaN(maxnuk,6*(maxnuk+1));
                end
            end
            C = NaN(maxnuk,6*(maxnuk+1));
            counter = 0;
            elCounter = 1;
            for ii = 1:d_p
                indices = indicesMat(ii,:);
                for jj = 1:size(plotAt,2)
                    if plotAt(ii,jj)
                        nuk1 = length(p_values{indices(2)});
                        nuk2 = length(p_values{indices(3)});
                        X_temp = zeros(nuk1, nuk2, 3);
                        dX_temp = zeros(nuk1, nuk2, 3, 3);
                        parfor j = 1:nuk1
    %                     for j = 1:nuk1
                            X_temp2 = zeros(1,nuk2, 3);
                            dX_temp2 = zeros(1,nuk2, 3, d_p);
                            Up_temp2 = cell(3,3);
                            for i_c = 1:d_p
                                for j_c = 1:d_p
                                    Up_temp2{i_c,j_c} = zeros(1,nuk2);
                                end
                            end
                            for k = 1:nuk2
                                xi = jj-1;
                                eta = p_values{indices(2)}(j);
                                zeta = p_values{indices(3)}(k);
                                XI = [xi,eta,zeta];
                                if plotJacobian || plotNormalVectors || plotParmDir
                                    [v,dvdxi,dvdeta,dvdzeta] = evaluateNURBS(nurbs, XI(indicesMat2(ii,:)),1);
                                    dX_temp2(1,k,:,1) = quiverScale*dvdxi/norm(dvdxi);
                                    dX_temp2(1,k,:,2) = quiverScale*dvdeta/norm(dvdeta);
                                    dX_temp2(1,k,:,3) = quiverScale*dvdzeta/norm(dvdzeta);
                                else
                                    v = evaluateNURBS(nurbs, XI(indicesMat2(ii,:)),1);
                                end
                                X_temp2(1,k,:) = v;
                            end
                            X_temp(j,:,:) = X_temp2;
                            dX_temp(j,:,:,:) = dX_temp2;
                        end
                        X(1:nuk1,counter+1:counter+nuk2,:) = X_temp;
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
                        if plotSolution
                            C(1:nuk1,counter+1:counter+nuk2) = reshape(colorFun(reshape(X_temp,nuk1*nuk2,3)),nuk1,nuk2);
                        end
                        
                        stepLen = res(indices(2:3))+1;
                        noLines = noUniqueKnots(indices(2:3));
                        
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

            maxC = max(max(C));
            minC = min(min(C));
            if plotSolution || plotJacobian
                surf(X(:,:,1),X(:,:,2),X(:,:,3),C,'EdgeColor','none','LineStyle','none','FaceAlpha',alphaValue, 'DisplayName',displayName)
                colorbar
            else
                surf(X(:,:,1),X(:,:,2),X(:,:,3), 'FaceColor', colorPatch,'EdgeColor','none','LineStyle','none','FaceAlpha',alphaValue, ...
                        'FaceLighting', faceLighting, 'DisplayName',displayName)
            end
            if plotParmDir
                for i = 1:d_p
                    quiver3(X(:,:,1),X(:,:,2),X(:,:,3),Up{1,i},Up{2,i},Up{3,i},'color',colorParmDirs{i},'AutoScale','off','DisplayName',[displayName, ' - parm dir ' num2str(i)])
                end
            end
            if plotElementEdges
                plotGridLines(vElementEdges,displayName)
            end
        elseif d_p == 2 && d == 3
            nuk1 = length(p_values{1});
            nuk2 = length(p_values{2});
            X = zeros(nuk1, nuk2, 3);
            dX = zeros(nuk1, nuk2, 3, d_p);
            normals = zeros(length(p_values{1}), length(p_values{2}),3);
            C = zeros(length(p_values{1}), length(p_values{2}));
            parfor i = 1:length(p_values{1})
                X_temp = zeros(length(p_values{2}),3);
                dX_temp = zeros(length(p_values{2}), 3, d_p);
                normals_temp = zeros(length(p_values{2}),3);
                C_temp = zeros(1,length(p_values{2}));
                xi = p_values{1}(i);
                for j = 1:length(p_values{2})
                    eta = p_values{2}(j);
                    [v,dvdxi,dvdeta] = evaluateNURBS(nurbs, [xi eta], 1);
                    dX_temp(j,:,1) = quiverScale*dvdxi/norm(dvdxi);
                    dX_temp(j,:,2) = quiverScale*dvdeta/norm(dvdeta);
                    normal = cross(dvdxi,dvdeta);
                    if plotSolution
                        ny = normal/norm(normal);
                        C_temp(j) = colorFun(v);
                    end
                    if plotNormalVectors
                        ny = quiverScale*normal/norm(normal);
                        normals_temp(j,:) = ny;
                    end
                    X_temp(j,:) = v;
                end
                X(i,:,:) = X_temp;
                dX(i,:,:,:) = dX_temp;
                normals(i,:,:) = normals_temp;
                C(i,:) = C_temp;
            end
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
            if plotNormalVectors
                quiver3(X(:,:,1),X(:,:,2),X(:,:,3),normals(:,:,1),normals(:,:,2),normals(:,:,3),'AutoScale','off','DisplayName',[displayName, ' - normal vectors'])
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
            if plotParmDir
                for i = 1:d_p
                    quiver3(X(:,:,1),X(:,:,2),X(:,:,3),dX(:,:,1,i),dX(:,:,2,i),dX(:,:,3,i),'color',colorParmDirs{i},'AutoScale','off','DisplayName',[displayName, ' - parm dir ' num2str(i)])
                end
            end
        elseif d_p == 2 && d == 2
            X = zeros(2*length(p_values{1})+2*length(p_values{2})-3,1);
            Y = zeros(2*length(p_values{1})+2*length(p_values{2})-3,1);
            counter = 1;
            % along eta = 0
            for i = 1:length(p_values{1})
                xi = p_values{1}(i);
                eta = 0;
                v = evaluateNURBS(nurbs, [xi eta]);
                X(counter) = v(1);
                Y(counter) = v(2);
                counter = counter + 1;
            end
            % along xi = 1
            for i = 2:length(p_values{2})
                xi = 1;
                eta = p_values{2}(i);
                v = evaluateNURBS(nurbs, [xi eta]);
                X(counter) = v(1);
                Y(counter) = v(2);
                counter = counter + 1;
            end
            % reverse order of arrays
            p_values{1} = p_values{1}(end:-1:1);
            p_values{2} = p_values{2}(end:-1:1);

            % along eta = 1
            for i = 2:length(p_values{1})
                xi = p_values{1}(i);
                eta = 1;
                v = evaluateNURBS(nurbs, [xi eta]);
                X(counter) = v(1);
                Y(counter) = v(2);
                counter = counter + 1;
            end
            % along xi = 0
            for i = 2:length(p_values{2})
                xi = 0;
                eta = p_values{2}(i);
                v = evaluateNURBS(nurbs, [xi eta]);
                X(counter) = v(1);
                Y(counter) = v(2);
                counter = counter + 1;
            end    
            % reverse back order of arrays
            p_values{1} = p_values{1}(end:-1:1);
            p_values{2} = p_values{2}(end:-1:1);
            fill(X,Y, colorPatch,'EdgeColor','none','LineStyle','none', 'DisplayName',displayName)
        elseif d_p == 1
            C = zeros(length(p_values{1}), d);

            parfor j = 1:length(p_values{1})
                xi = p_values{1}(j);
                C(j,:) = evaluateNURBS(nurbs, xi);
            end
            switch d
                case 1
                    plot(C,zeros(size(C)), 'color', colorPatch, 'DisplayName',displayName);  
                case 2
                    plot(C(:,1), C(:,2), 'color', colorPatch, 'DisplayName',displayName);  
                case 3
                    plot3(C(:,1), C(:,2), C(:,3), 'color', colorPatch, 'DisplayName',displayName); 
            end  
        end
    end
    if plotElementEdges && d_p == 2 && d == 2
        v = [];
        for i = 1:d_p
            indices = indicesMat(i,:);
            nuk = length(uniqueKnots{indices(2)}); % number of unique knots
            if d_p == 3
                nuk2 = length(uniqueKnots{indices(3)}); % number of unique knots
            else
                nuk2 = 1;
            end
            nok = length(p_values{indices(1)}); % number of knots
            temp = zeros(d,nuk*nuk2*(nok+1));
            counter = 0;
            for j = 1:nuk
                eta = uniqueKnots{indices(2)}(j);
                if d_p == 3
                    for k = 1:nuk2
                        zeta = uniqueKnots{indices(3)}(k);
                        temp2 = zeros(d,nok);
                        parfor ii = 1:nok
                            xi = p_values{indices(1)}(ii);
                            XI = [xi eta zeta];
                            temp2(:,ii) = evaluateNURBS(nurbs, XI(indicesMat2(i,:)));
                        end
                        temp(:,counter+1:counter+nok+1) = [temp2, NaN(d,1)];
                        counter = counter + nok+1;
                    end
                else
                    temp2 = zeros(d,nok);
                    parfor ii = 1:nok
                        xi = p_values{indices(1)}(ii);
                        XI = [xi eta];
                        temp2(:,ii) = evaluateNURBS(nurbs, XI(indicesMat2(i,:)));
                    end
                    temp(:,counter+1:counter+nok+1) = [temp2, NaN(d,1)];
                    counter = counter + nok+1;
                end
            end
            v = [v, temp];
            if size(v,1) > 2
                plot3(v(1,:),v(2,:),v(3,:),'color',lineColor,'LineWidth',lineWidth,'DisplayName',[displayName, ' - element edges'])
            else
                plot(v(1,:),v(2,:),'color',lineColor,'LineWidth',lineWidth,'DisplayName',[displayName, ' - element edges'])
            end
        end
    end
    if plotControlPolygon
        if d_p == 1
            v = coeffs(1:d,:);
        else
            v = [];
            for i = 1:d_p
                indices = 1:d_p+1;
                indices([i+1,2]) = [2,i+1];
                temp = permute(coeffs,indices);
                prd = prod(dimensions(indices(3:end))); % product of remaining dimensions
                temp = cat(2,reshape(temp,dimensions(1),dimensions(i+1),prd),NaN(d,1,prd));  % adding NaN to make the data discontinous
                temp = reshape(temp,dimensions(1),(dimensions(i+1)+1)*prd);
                v = [v, temp, NaN(dimensions(1),1)];
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
    view(options.view)
end
axis equal
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
