function [h,maxC] = plotNURBS(nurbs, newOptions)

%% Interpret input arguments

% set default values
options = struct('resolution',[10,10,10],...
                 'plotElementEdges',true,...
                 'plotControlPolygon', true, ...
                 'plotNormalVectors', false, ...
                 'color',1.5*[44 77 32]/255,...
                 'alphaValue',1,...
                 'plotAt',ones(3,2),...
                 'colorFun', NaN,...
                 'lineColor', 'black', ...
                 'colorControlPolygon', 'red', ...
                 'markerEdgeColor', 'black', ...
                 'displayName', '', ...
                 'markerColor', 'black', ...
                 'samplingDistance', NaN, ...
                 'LineWidth',0.5, ...
                 'elementBasedSamples',false);
options = updateOptions(options,newOptions);


sd = options.samplingDistance;
resolution = options.resolution;
plotElementEdges = options.plotElementEdges;
displayName = options.displayName;
plotControlPolygon = options.plotControlPolygon;
plotNormalVectors = options.plotNormalVectors;
markerEdgeColor = options.markerEdgeColor;
color = options.color;
alphaValue = options.alphaValue;
plotAt = options.plotAt;
colorFun = options.colorFun;
lineWidth = options.LineWidth;
lineColor = options.lineColor;
colorControlPolygon = options.colorControlPolygon;
markerColor = options.markerColor;
elementBasedSamples = options.elementBasedSamples;

plotSolution = isa(colorFun, 'function_handle');


coeffs = nurbs.coeffs;
d_p = numel(nurbs.knots);
dimensions = size(coeffs);
d = dimensions(1)-1;
p_values = cell(1,d_p);
uniqueKnots = cell(1,d_p);
for i = 1:d_p
    % To reconstruct any interpolatory points, any p repeated knot should be
    % included in the abscissa values to plot
    uniqueKnots{i} = unique(nurbs.knots{i});
    p_values{i} = insertUniform3(nurbs.knots{i}, resolution(i), nurbs.degree(i));
end
hold on
switch nurbs.type
    case '3Dvolume'
        indices = 1:d_p;
        nuk1 = length(p_values{indices(1)});
        nuk2 = length(p_values{indices(2)});
        nuk3 = length(p_values{indices(3)});
        maxnuk = max([nuk1,nuk2,nuk3]);
        X = NaN(maxnuk,6*(maxnuk+1));
        Y = NaN(maxnuk,6*(maxnuk+1));
        Z = NaN(maxnuk,6*(maxnuk+1));
        C = NaN(maxnuk,6*(maxnuk+1));
        counter = 0;
        for ii = 1:d_p
            indices = indices([2:end,1]);
            for jj = 1:size(plotAt,2)
                %% Plot surface where xi = 0
                if plotAt(ii,jj)
                    nuk1 = length(p_values{indices(2)});
                    nuk2 = length(p_values{indices(3)});
                    X_temp = zeros(nuk1, nuk2);
                    Y_temp = zeros(nuk1, nuk2);
                    Z_temp = zeros(nuk1, nuk2);
                    C_temp = zeros(nuk1, nuk2);
                    parfor j = 1:nuk1
%                     for j = 1:nuk1
                        X_temp2 = zeros(1,nuk2);
                        Y_temp2 = zeros(1,nuk2);
                        Z_temp2 = zeros(1,nuk2);
                        C_temp2 = zeros(1,nuk2);
                        for k = 1:nuk2
                            xi = jj-1;
                            eta = p_values{indices(2)}(j);
                            zeta = p_values{indices(3)}(k);
                            XI = [xi,eta,zeta];
                            v = evaluateNURBS(nurbs, XI(indices));
                            X_temp2(k) = v(1);
                            Y_temp2(k) = v(2);
                            Z_temp2(k) = v(3);
                            if plotSolution
                                C_temp2(k) = colorFun(v);
                            end
                        end
                        X_temp(j,:) = X_temp2;
                        Y_temp(j,:) = Y_temp2;
                        Z_temp(j,:) = Z_temp2;
                        C_temp(j,:) = C_temp2;
                    end
                    X(1:nuk1,counter+1:counter+nuk2) = X_temp;
                    Y(1:nuk1,counter+1:counter+nuk2) = Y_temp;
                    Z(1:nuk1,counter+1:counter+nuk2) = Z_temp;
                    C(1:nuk1,counter+1:counter+nuk2) = C_temp;
                    counter = counter + nuk2+1;
                end
            end
        end
        if plotSolution
            if plotElementEdges
                surf(X,Y,Z, C, 'EdgeColor','none','LineStyle','none', 'DisplayName',displayName)
            else
                surf(X,Y,Z, C)
            end
        else
            if true
                if alphaValue < 1
                    surf(X,Y,Z, 'FaceColor', color,'EdgeColor','none','LineStyle','none','FaceAlpha',alphaValue, ...
                                'FaceLighting', 'none', 'DisplayName',displayName)
                else
                    surf(X,Y,Z, 'FaceColor', color,'EdgeColor','none','LineStyle','none','FaceAlpha',alphaValue, ...
                                'DisplayName',displayName)
                end
            else
                surf(X,Y,Z, 'FaceColor', color, 'DisplayName',displayName)
            end
        end
    case '3Dsurface'
        if 0 %plotSolution
            unif = linspace(-1, 1, resolution(1));
            p_values{1} = insertUniform3(nurbs.knots{1}, resolution(1), nurbs.degree(1));
            xi_x = varCol.xi_x;
            if xi_x < 0.2
                coeff = xi_x;
            elseif xi_x > 0.8
                coeff = 1-xi_x;
            else
                coeff = 0.2;
            end
            p_values{1} = sort([p_values{1}', xi_x+coeff*sign(unif).*abs(unif).^1.3])';

            unif = linspace(-1, 1, resolution(2));
            p_values{2} = insertUniform3(nurbs.knots{2}, resolution(2), nurbs.degree(2));
            eta_x = varCol.eta_x;
            if eta_x < 0.2
                coeff = eta_x;
            elseif eta_x > 0.8
                coeff = 1-eta_x;
            else
                coeff = 0.2;
            end
            p_values{2} = sort([p_values{2}', eta_x+coeff*sign(unif).*abs(unif).^1.3])';
        else
            if elementBasedSamples
                p_values{1} = insertUniform(nurbs.knots{1}, resolution(1));
                p_values{2} = insertUniform(nurbs.knots{2}, resolution(2));
            else
                p_values{1} = insertUniform3(nurbs.knots{1}, resolution(1), nurbs.degree(1));
                p_values{2} = insertUniform3(nurbs.knots{2}, resolution(2), nurbs.degree(2));
            end
        end
    %     %% Plot surface where zeta = 1
        X = zeros(length(p_values{1}), length(p_values{2}));
        Y = zeros(length(p_values{1}), length(p_values{2}));
        Z = zeros(length(p_values{1}), length(p_values{2}));
        U = zeros(length(p_values{1}), length(p_values{2}));
        V = zeros(length(p_values{1}), length(p_values{2}));
        W = zeros(length(p_values{1}), length(p_values{2}));
        C = zeros(length(p_values{1}), length(p_values{2}));
        if 0 %plotSolution
            p_xi = nurbs.degree(1);
            p_eta = nurbs.degree(2);
            Xi = nurbs.knots{1};
            Eta = nurbs.knots{2};
            weights = varCol.weights;
            controlPts = varCol.controlPts;
            colorFun = varCol.colorFun;
            x = varCol.x;
            d_vec = varCol.d_vec;
        else
            p_xi = NaN;
            p_eta = NaN;
            Xi = NaN;
            Eta = NaN;
            weights = NaN;
            controlPts = NaN;
            x = NaN;
            d_vec = NaN;
        end
        if isnan(sd)
%             for i = 1:length(p_values{1})
            parfor i = 1:length(p_values{1})
                X_temp = zeros(1,length(p_values{2}));
                Y_temp = zeros(1,length(p_values{2}));
                Z_temp = zeros(1,length(p_values{2}));
                U_temp = zeros(1,length(p_values{2}));
                V_temp = zeros(1,length(p_values{2}));
                W_temp = zeros(1,length(p_values{2}));
                C_temp = zeros(1,length(p_values{2}));
                xi = p_values{1}(i);
                for j = 1:length(p_values{2})
                    eta = p_values{2}(j);
                    [y,dydxi,dydeta] = evaluateNURBS_2ndDeriv(nurbs, [xi eta]);
                    normal = cross(dydxi,dydeta);
                    if plotSolution
        %                 [u, v] = numericalSolEval_final_surf(xi, eta, p_xi, p_eta, Xi, Eta, weights, controlPts, U);
    %                     y = evaluateNURBS(nurbs, [xi eta]);
    %                     xmy = x - y.';
    %                     nx = x'/norm(x');
    %                     ny = y/norm(y);
    %                     r = norm(xmy);
        %                 C_temp(j) = real(colorFun(xmy,r,nx,ny));
        %                 C_temp(j) = log10(abs(norm(-xmy/r + d_vec')));
    %                     C_temp(j) = colorFun(y);
    %                     C_temp(j) = abs(real(colorFun(y.')));
                        [y,dydxi,dydeta] = evaluateNURBS_2ndDeriv(nurbs, [xi eta]);
                        ny = normal/norm(normal);
                        dotProd = dot([0,-1,0],ny);
                        k = 2*pi*1000/1500;
                        if false
                            C_temp2 = 0;
                            if dotProd > 0
                                C_temp2 = real(dotProd*exp(-2*1i*k*dot([0,-1,0],y)));
                            end
                            C_temp(j) = C_temp2;
                        else
                            C_temp(j) = colorFun(y);
                        end
    %                     C_temp(j) = colorFun(y.',ny);
    %                     y = y + ny*real(C_temp(j))*2;
        %                 C_temp(j)
        %                 C_temp(j) = colorFun(y,xi,eta);
        %                 C_temp(j) = log10(abs(norm(y)-1));
        %                 C_temp(j) = abs(norm(y)-1);
                    end
                    X_temp(j) = y(1);
                    Y_temp(j) = y(2);
                    Z_temp(j) = y(3);
                    U_temp(j) = normal(1);
                    V_temp(j) = normal(2);
                    W_temp(j) = normal(3);
                end
                X(i,:) = X_temp;
                Y(i,:) = Y_temp;
                Z(i,:) = Z_temp;
                U(i,:) = U_temp;
                V(i,:) = V_temp;
                W(i,:) = W_temp;
                C(i,:) = C_temp;
            end
            if plotSolution
                if false
                    noPts = length(p_values{1})*length(p_values{2});
                    C = reshape(colorFun([reshape(X,noPts,1),reshape(Y,noPts,1),reshape(Z,noPts,1)],...
                                         [reshape(X,noPts,1),reshape(Y,noPts,1),reshape(Z,noPts,1)].'), ...
                                          length(p_values{1}), length(p_values{2}));
        %             C = log10(abs(C-analytic_v)/max(max(abs(C))));
                end
                if plotElementEdges
                    surf(X,Y,Z, C, 'EdgeColor','none','LineStyle','none','DisplayName',displayName)
                else
                    surf(X,Y,Z, C,'DisplayName',displayName)
                end
                maxC = max(max(C));
                colorbar
            else
                if plotNormalVectors
                    quiver3(X,Y,Z,U,V,W,'DisplayName',[displayName, ' - normal vectors'])
                end
                if plotElementEdges
                    if alphaValue < 1
                        surf(X,Y,Z, 'FaceColor', color,'EdgeColor','none','LineStyle','none','FaceAlpha',alphaValue, ...
                                    'FaceLighting', 'none','DisplayName',displayName)
                    else
                        surf(X,Y,Z, 'FaceColor', color,'EdgeColor','none','LineStyle','none','FaceAlpha',alphaValue,...
                                    'DisplayName',displayName)
                    end
                else
                    surf(X,Y,Z, 'FaceColor', color,'DisplayName',displayName)
                end
            end
        else
            varCol.dimension = 1;
            varCol.nurbs = nurbs;
            varCol = convertNURBS(varCol);  
            varCol = generateIGA2DMesh_new(varCol);
            noElems = varCol.patches{1}.noElems;
            index = varCol.patches{1}.index;
            elRangeXi = varCol.patches{1}.elRange{1};
            elRangeEta = varCol.patches{1}.elRange{2};
            maxC = zeros(noElems,1);
            X = cell(noElems,1);
            Y = cell(noElems,1);
            Z = cell(noElems,1);
            C = cell(noElems,1);
            parfor e = 1:noElems
%             for e = 1:noElems
                idXi = index(e,1);
                idEta = index(e,2);

                Xi_e = elRangeXi(idXi,:);
                Eta_e = elRangeEta(idEta,:);
                [~, dydxi, dydeta] = evaluateNURBS_2ndDeriv(nurbs, [Xi_e(1) Eta_e(1); Xi_e(2) Eta_e(1); Xi_e(2) Eta_e(2); Xi_e(1) Eta_e(2)]);
                d12 = mean([norm(dydxi(1,:)),norm(dydxi(2,:))])*(Xi_e(2)-Xi_e(1));
                d43 = mean([norm(dydxi(4,:)),norm(dydxi(3,:))])*(Xi_e(2)-Xi_e(1));
                dxi = max(d12,d43);
                d23 = mean([norm(dydeta(2,:)),norm(dydeta(3,:))])*(Eta_e(2)-Eta_e(1));
                d14 = mean([norm(dydeta(1,:)),norm(dydeta(4,:))])*(Eta_e(2)-Eta_e(1));
                deta = max(d23,d14);


                noXiValues  = max([round(dxi/sd)+1, 2]);
                noEtaValues = max([round(deta/sd)+1, 2]);
                if noXiValues*noEtaValues < resolution(1)*resolution(2)
                    noXiValues = max([round(resolution(1)*sqrt(dxi/deta)),2]);
                    noEtaValues = max([round(resolution(2)*sqrt(deta/dxi)),2]);
                end
                Xi_values2 = linspace(Xi_e(1),Xi_e(2),noXiValues).';
                Eta_values2 = linspace(Eta_e(1),Eta_e(2),noEtaValues).';
                xi = copyVector(Xi_values2,noEtaValues,1);
                eta = copyVector(Eta_values2,noXiValues,2);

                v = evaluateNURBS_2ndDeriv(nurbs, [xi eta]);

                X{e} = reshape(v(:,1),noXiValues,noEtaValues);
                Y{e} = reshape(v(:,2),noXiValues,noEtaValues);
                Z{e} = reshape(v(:,3),noXiValues,noEtaValues); 
                if plotSolution
                    C{e} = reshape(colorFun(v),noXiValues,noEtaValues);
        %             C{e} = zeros(noXiValues,noEtaValues);
                    maxC(e) = max(max(C{e}));
                else
                    maxC(e) = NaN;
                end
            end
            for e = 1:noElems
                if plotSolution
                    if plotElementEdges
                        surf(X{e},Y{e},Z{e}, C{e}, 'EdgeColor','none','LineStyle','none','DisplayName',displayName)
                    else
                        surf(X{e},Y{e},Z{e}, C{e},'DisplayName',displayName)
                    end
                else
                    if plotElementEdges
                        if alphaValue < 1
                            surf(X{e},Y{e},Z{e}, 'FaceColor', color,'EdgeColor','none','LineStyle','none', ...
                                'FaceAlpha',alphaValue, 'FaceLighting', 'none','DisplayName',displayName)
                        else
                            surf(X{e},Y{e},Z{e}, 'FaceColor', color,'EdgeColor','none','LineStyle','none', ...
                                'FaceAlpha',alphaValue,'DisplayName',displayName)
                        end
                    else
                        surf(X{e},Y{e},Z{e})
                    end
                end
            end
            if plotSolution
                maxC = max(maxC);
                colorbar
            end
        end
    case '2Dsurface'
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
        if plotElementEdges
            fill(X,Y, color,'EdgeColor','none','LineStyle','none')
        else
            fill(X,Y, color)
        end
    case '3Dcurve'
        p_values{1} = insertUniform(nurbs.knots{1}, resolution);

        C = zeros(length(p_values{1}), 3);

        for j = 1:length(p_values{1})
            xi = p_values{1}(j);
            C(j,:) = evaluateNURBS(nurbs, xi);
        end

        plot3(C(:,1), C(:,2), C(:,3), 'color', color); 
    case '2Dcurve'
        p_values{1} = insertUniform(nurbs.knots{1}, resolution);

        C = zeros(length(p_values{1}), 2);

        for j = 1:length(p_values{1})
            xi = p_values{1}(j);
            C(j,:) = evaluateNURBS(nurbs, xi);
        end

        plot(C(:,1), C(:,2), 'color', color);  
    case '1Dnurbs'
        p_values{1} = insertUniform(nurbs.knots{1}, resolution);

        C = zeros(length(p_values{1}), 1);

        for j = 1:length(p_values{1})
            xi = p_values{1}(j);
            C(j) = evaluateNURBS(nurbs, xi);
        end

        plot(p_values{1}, C, 'color', color);    
end
if plotControlPolygon
    v = [];
    for i = 1:d_p
        indices = 1:d_p+1;
        indices([i+1,2]) = [2,i+1];
        temp = permute(coeffs,indices);
        prd = prod(dimensions(indices(3:end))); % product of remaining dimensions
        temp = reshape(temp,dimensions(1),dimensions(i+1),prd);
        temp(:,end+1,:) = NaN; % make data discontinous
        temp = reshape(temp,dimensions(1),(dimensions(i+1)+1)*prd);
        v = [v, temp, NaN(dimensions(1),1)];
    end
    if size(v,1) > 2
        plot3(v(1,:),v(2,:),v(3,:),'o-','color',colorControlPolygon,'MarkerFaceColor', markerColor, ...
                                        'MarkerEdgeColor', markerEdgeColor,'DisplayName',[displayName, ' - control polygon'])
    else
        plot(v(1,:),v(2,:),'o-','color',colorControlPolygon,'MarkerFaceColor', markerColor, ...
                                'MarkerEdgeColor', markerEdgeColor,'DisplayName',[displayName, ' - control polygon'])
    end
end


if plotElementEdges   
    indices = 1:d_p;
    v = [];
    for i = 1:d_p
        indices = indices([2:end,1]);
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
                        XI = [xi eta, zeta];
                        temp2(:,ii) = evaluateNURBS(nurbs, XI(indices));
                    end
                    temp(:,counter+1:counter+nok+1) = [temp2, NaN(d,1)];
                    counter = counter + nok+1;
                end
            else
                temp2 = zeros(d,nok);
                parfor ii = 1:nok
                    xi = p_values{indices(1)}(ii);
                    XI = [xi eta];
                    temp2(:,ii) = evaluateNURBS(nurbs, XI(indices));
                end
                temp(:,counter+1:counter+nok+1) = [temp2, NaN(d,1)];
                counter = counter + nok+1;
            end
        end
        v = [v, temp];
    end
    if size(v,1) > 2
        plot3(v(1,:),v(2,:),v(3,:),'color',lineColor,'LineWidth',lineWidth,'DisplayName',[displayName, ' - element edges'])
    else
        plot(v(1,:),v(2,:),'color',lineColor,'LineWidth',lineWidth,'DisplayName',[displayName, ' - element edges'])
    end
end
h = gcf;


