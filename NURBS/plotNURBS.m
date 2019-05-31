function [h,maxC] = plotNURBS(nurbs, newOptions)

%% Interpret input arguments

% set default values
options = struct('resolution',[10,10,10],...
                 'plotElementEdges',true,...
                 'color',1.5*[44 77 32]/255,...
                 'alphaValue',1,...
                 'plotAt',ones(3,2),...
                 'colorFun', NaN,...
                 'lineColor', 'black', ...
                 'samplingDistance', NaN, ...
                 'LineWidth',0.5, ...
                 'elementBasedSamples',false);

% read the acceptable names
optionNames = fieldnames(options);

if nargin == 2
    % count arguments
    nArgs = length(newOptions);
    if round(nArgs/2) ~= nArgs/2
        error('Must have propertyName/propertyValue pairs')
    end

    for pair = reshape(newOptions,2,[]) %# pair is {propName;propValue}
        inpName = pair{1}; %# make case insensitive

        if any(strcmp(inpName,optionNames))
            options.(inpName) = pair{2};
        else
            error('%s is not a recognized parameter name',inpName)
        end
    end
end


d = options.samplingDistance;
resolution = options.resolution;
plotElementEdges = options.plotElementEdges;
color = options.color;
alphaValue = options.alphaValue;
plotAt = options.plotAt;
colorFun = options.colorFun;
lineWidth = options.LineWidth;
lineColor = options.lineColor;
elementBasedSamples = options.elementBasedSamples;

plotSolution = isa(colorFun, 'function_handle');

if strcmp(nurbs.type, '3Dvolume')
    unique_xiKnots = unique(nurbs.knots{1});
    unique_etaKnots = unique(nurbs.knots{2});
    unique_zetaKnots = unique(nurbs.knots{3});
    
    % To reconstruct any interpolatory points, any p repeated knot should be
    % included in the abscissa values to plot
    Xi_values = insertUniform3(nurbs.knots{1}, resolution(1), nurbs.degree(1));
    Eta_values = insertUniform3(nurbs.knots{2}, resolution(2), nurbs.degree(2));
    Zeta_values = insertUniform3(nurbs.knots{3}, resolution(3), nurbs.degree(3));

    %% Plot surface where xi = 0
    if plotAt(1,1)
        X = zeros(length(Eta_values), length(Zeta_values));
        Y = zeros(length(Eta_values), length(Zeta_values));
        Z = zeros(length(Eta_values), length(Zeta_values));
        C = zeros(length(Eta_values), length(Zeta_values));
        parfor j = 1:length(Eta_values)
%         for j = 1:length(Eta_values)
            X_temp = zeros(1,length(Zeta_values));
            Y_temp = zeros(1,length(Zeta_values));
            Z_temp = zeros(1,length(Zeta_values));
            C_temp = zeros(1,length(Zeta_values));
            for k = 1:length(Zeta_values)
                xi = 0;
                eta = Eta_values(j);
                zeta = Zeta_values(k);
                v = evaluateNURBS(nurbs, [xi eta zeta]);
                X_temp(k) = v(1);
                Y_temp(k) = v(2);
                Z_temp(k) = v(3);
                if plotSolution
                    C_temp(k) = colorFun(v);
                end
            end
            X(j,:) = X_temp;
            Y(j,:) = Y_temp;
            Z(j,:) = Z_temp;
            C(j,:) = C_temp;
        end
        if plotSolution
            if plotElementEdges
                surf(X,Y,Z, C, 'EdgeColor','none','LineStyle','none')
            else
                surf(X,Y,Z, C)
            end
        else
            if true
                if alphaValue < 1
                    surf(X,Y,Z, 'FaceColor', color,'EdgeColor','none','LineStyle','none','FaceAlpha',alphaValue, 'FaceLighting', 'none')
                else
                    surf(X,Y,Z, 'FaceColor', color,'EdgeColor','none','LineStyle','none','FaceAlpha',alphaValue)
                end
            else
                surf(X,Y,Z, 'FaceColor', color)
            end
        end
        hold on
    end


    % Plot surface where xi = 1
    if plotAt(1,2)
        X = zeros(length(Eta_values), length(Zeta_values));
        Y = zeros(length(Eta_values), length(Zeta_values));
        Z = zeros(length(Eta_values), length(Zeta_values));
        C = zeros(length(Eta_values), length(Zeta_values));
        parfor j = 1:length(Eta_values)
            X_temp = zeros(1,length(Zeta_values));
            Y_temp = zeros(1,length(Zeta_values));
            Z_temp = zeros(1,length(Zeta_values));
            C_temp = zeros(1,length(Zeta_values));
            for k = 1:length(Zeta_values)
                xi = 1;
                eta = Eta_values(j);
                zeta = Zeta_values(k);
                v = evaluateNURBS(nurbs, [xi eta zeta]);
                X_temp(k) = v(1);
                Y_temp(k) = v(2);
                Z_temp(k) = v(3);
                if plotSolution
                    C_temp(k) = colorFun(v);
                end
            end
            X(j,:) = X_temp;
            Y(j,:) = Y_temp;
            Z(j,:) = Z_temp;
            C(j,:) = C_temp;
        end
        if plotSolution
            if plotElementEdges
                surf(X,Y,Z, C, 'EdgeColor','none','LineStyle','none')
            else
                surf(X,Y,Z, C)
            end
        else
            if true
                if alphaValue < 1
                    surf(X,Y,Z, 'FaceColor', color,'EdgeColor','none','LineStyle','none','FaceAlpha',alphaValue, 'FaceLighting', 'none')
                else
                    surf(X,Y,Z, 'FaceColor', color,'EdgeColor','none','LineStyle','none','FaceAlpha',alphaValue)
                end
            else
                surf(X,Y,Z, 'FaceColor', color)
            end
        end
        hold on
    end
    

    % Plot surface where eta = 0
    if plotAt(2,1)
        X = zeros(length(Xi_values), length(Zeta_values));
        Y = zeros(length(Xi_values), length(Zeta_values));
        Z = zeros(length(Xi_values), length(Zeta_values));
        C = zeros(length(Xi_values), length(Zeta_values));
        parfor i = 1:length(Xi_values)
            X_temp = zeros(1,length(Zeta_values));
            Y_temp = zeros(1,length(Zeta_values));
            Z_temp = zeros(1,length(Zeta_values));
            C_temp = zeros(1,length(Zeta_values));
            for k = 1:length(Zeta_values)
                xi = Xi_values(i);
                eta = 0;
                zeta = Zeta_values(k);
                v = evaluateNURBS(nurbs, [xi eta zeta]);
                X_temp(k) = v(1);
                Y_temp(k) = v(2);
                Z_temp(k) = v(3);
                if plotSolution
                    C_temp(k) = colorFun(v);
                end
            end
            X(i,:) = X_temp;
            Y(i,:) = Y_temp;
            Z(i,:) = Z_temp;
            C(i,:) = C_temp;
        end

        if plotSolution
            if plotElementEdges
                surf(X,Y,Z, C, 'EdgeColor','none','LineStyle','none')
            else
                surf(X,Y,Z, C)
            end
        else
            if true
                if alphaValue < 1
                    surf(X,Y,Z, 'FaceColor', color,'EdgeColor','none','LineStyle','none','FaceAlpha',alphaValue, 'FaceLighting', 'none')
                else
                    surf(X,Y,Z, 'FaceColor', color,'EdgeColor','none','LineStyle','none','FaceAlpha',alphaValue)
                end
            else
                surf(X,Y,Z, 'FaceColor', color)
            end
        end
        hold on
    end
    

    if plotAt(2,2)
        % Plot surface where eta = 1
        X = zeros(length(Xi_values), length(Zeta_values));
        Y = zeros(length(Xi_values), length(Zeta_values));
        Z = zeros(length(Xi_values), length(Zeta_values));
        C = zeros(length(Xi_values), length(Zeta_values));
        parfor i = 1:length(Xi_values)
            X_temp = zeros(1,length(Zeta_values));
            Y_temp = zeros(1,length(Zeta_values));
            Z_temp = zeros(1,length(Zeta_values));
            C_temp = zeros(1,length(Zeta_values));
            for k = 1:length(Zeta_values)
                xi = Xi_values(i);
                eta = 1;
                zeta = Zeta_values(k);
                v = evaluateNURBS(nurbs, [xi eta zeta]);
                X_temp(k) = v(1);
                Y_temp(k) = v(2);
                Z_temp(k) = v(3);
                if plotSolution
                    C_temp(k) = colorFun(v);
                end
            end
            X(i,:) = X_temp;
            Y(i,:) = Y_temp;
            Z(i,:) = Z_temp;
            C(i,:) = C_temp;
        end
        if plotSolution
            if plotElementEdges
                surf(X,Y,Z, C, 'EdgeColor','none','LineStyle','none')
            else
                surf(X,Y,Z, C)
            end
        else
            if true
                if alphaValue < 1
                    surf(X,Y,Z, 'FaceColor', color,'EdgeColor','none','LineStyle','none','FaceAlpha',alphaValue, 'FaceLighting', 'none')
                else
                    surf(X,Y,Z, 'FaceColor', color,'EdgeColor','none','LineStyle','none','FaceAlpha',alphaValue)
                end
            else
                surf(X,Y,Z, 'FaceColor', color)
            end
        end
        hold on
    end
    

    if plotAt(3,1)
        %% Plot surface where zeta = 0
        X = zeros(length(Xi_values), length(Eta_values));
        Y = zeros(length(Xi_values), length(Eta_values));
        Z = zeros(length(Xi_values), length(Eta_values));
        C = zeros(length(Xi_values), length(Eta_values));
        parfor i = 1:length(Xi_values)
            X_temp = zeros(1,length(Eta_values));
            Y_temp = zeros(1,length(Eta_values));
            Z_temp = zeros(1,length(Eta_values));
            C_temp = zeros(1,length(Eta_values));
            for j = 1:length(Eta_values)
                xi = Xi_values(i);
                eta = Eta_values(j);
                zeta = 0;
                v = evaluateNURBS(nurbs, [xi eta zeta]);
                X_temp(j) = v(1);
                Y_temp(j) = v(2);
                Z_temp(j) = v(3);
                if plotSolution
                    C_temp(j) = colorFun(v);
                end
            end
            X(i,:) = X_temp;
            Y(i,:) = Y_temp;
            Z(i,:) = Z_temp;
            C(i,:) = C_temp;
        end
        if plotSolution
            if plotElementEdges
                surf(X,Y,Z, C, 'EdgeColor','none','LineStyle','none')
            else
                surf(X,Y,Z, C)
            end
        else
            if true
                if alphaValue < 1
                    surf(X,Y,Z, 'FaceColor', color,'EdgeColor','none','LineStyle','none','FaceAlpha',alphaValue, 'FaceLighting', 'none')
                else
                    surf(X,Y,Z, 'FaceColor', color,'EdgeColor','none','LineStyle','none','FaceAlpha',alphaValue)
                end
            else
                surf(X,Y,Z, 'FaceColor', color)
            end
        end
        hold on
    end

    if plotAt(3,2)
        %% Plot surface where zeta = 1
        X = zeros(length(Xi_values), length(Eta_values));
        Y = zeros(length(Xi_values), length(Eta_values));
        Z = zeros(length(Xi_values), length(Eta_values));
        C = zeros(length(Xi_values), length(Eta_values));
        parfor i = 1:length(Xi_values)
            X_temp = zeros(1,length(Eta_values));
            Y_temp = zeros(1,length(Eta_values));
            Z_temp = zeros(1,length(Eta_values));
            C_temp = zeros(1,length(Eta_values));
            for j = 1:length(Eta_values)
                xi = Xi_values(i);
                eta = Eta_values(j);
                zeta = 1;
                v = evaluateNURBS(nurbs, [xi eta zeta]);
                X_temp(j) = v(1);
                Y_temp(j) = v(2);
                Z_temp(j) = v(3);
                if plotSolution
                    C_temp(j) = colorFun(v);
                end
            end
            X(i,:) = X_temp;
            Y(i,:) = Y_temp;
            Z(i,:) = Z_temp;
            C(i,:) = C_temp;
        end
        if plotSolution
            if plotElementEdges
                surf(X,Y,Z, C, 'EdgeColor','none','LineStyle','none')
            else
                surf(X,Y,Z, C)
            end
            colorbar
        else
            if true
                if alphaValue < 1
                    surf(X,Y,Z, 'FaceColor', color,'EdgeColor','none','LineStyle','none','FaceAlpha',alphaValue, 'FaceLighting', 'none')
                else
                    surf(X,Y,Z, 'FaceColor', color,'EdgeColor','none','LineStyle','none','FaceAlpha',alphaValue)
                end
            else
                surf(X,Y,Z, 'FaceColor', color)
            end
        end
        hold on
    end

    if plotElementEdges       
        x = zeros(length(Xi_values),1);
        y = zeros(length(Xi_values),1);
        z = zeros(length(Xi_values),1);
     
        for j = 1:length(unique_etaKnots)
            eta = unique_etaKnots(j);
            for k = 1:length(unique_zetaKnots)
                zeta = unique_zetaKnots(k);
                for i = 1:length(Xi_values)
                    xi = Xi_values(i);
                    v = evaluateNURBS(nurbs, [xi eta zeta]);
                    x(i) = v(1);
                    y(i) = v(2);
                    z(i) = v(3);
                end
                plot3(x,y,z,'color',lineColor)
            end
        end
        
        x = zeros(length(Eta_values),1);
        y = zeros(length(Eta_values),1);
        z = zeros(length(Eta_values),1);
        for i = 1:length(unique_xiKnots)
            xi = unique_xiKnots(i);
            for k = 1:length(unique_zetaKnots)
                zeta = unique_zetaKnots(k);
                for j = 1:length(Eta_values)
                    eta = Eta_values(j);
                    v = evaluateNURBS(nurbs, [xi eta zeta]);
                    x(j) = v(1);
                    y(j) = v(2);
                    z(j) = v(3);
                end
                plot3(x,y,z,'color',lineColor)
            end
        end

        x = zeros(length(Zeta_values),1);
        y = zeros(length(Zeta_values),1);
        z = zeros(length(Zeta_values),1);
        for i = 1:length(unique_xiKnots)
            xi = unique_xiKnots(i);
            for j = 1:length(unique_etaKnots)
                eta = unique_etaKnots(j);
                for k = 1:length(Zeta_values)
                    zeta = Zeta_values(k);
                    v = evaluateNURBS(nurbs, [xi eta zeta]);
                    x(k) = v(1);
                    y(k) = v(2);
                    z(k) = v(3);
                end
                plot3(x,y,z,'color',lineColor)
            end
        end
    end
% 
%     daspect([1 1 1]);
%     view(3)
%     camlight(110,30);
%     lighting phong;
%     set(gca, 'Color', 'none');
%     camlight('left')
    h = gcf;
elseif strcmp(nurbs.type, '3Dsurface')
    unique_xiKnots = unique(nurbs.knots{1});
    unique_etaKnots = unique(nurbs.knots{2});
    
    % To reconstruct any interpolatory points, any p repeated knot should be
    % included in the abscissa values to plot
    if 0 %plotSolution
        unif = linspace(-1, 1, resolution(1));
        Xi_values = insertUniform3(nurbs.knots{1}, resolution(1), nurbs.degree(1));
        xi_x = varCol.xi_x;
        if xi_x < 0.2
            coeff = xi_x;
        elseif xi_x > 0.8
            coeff = 1-xi_x;
        else
            coeff = 0.2;
        end
        Xi_values = sort([Xi_values', xi_x+coeff*sign(unif).*abs(unif).^1.3])';
        
        unif = linspace(-1, 1, resolution(2));
        Eta_values = insertUniform3(nurbs.knots{2}, resolution(2), nurbs.degree(2));
        eta_x = varCol.eta_x;
        if eta_x < 0.2
            coeff = eta_x;
        elseif eta_x > 0.8
            coeff = 1-eta_x;
        else
            coeff = 0.2;
        end
        Eta_values = sort([Eta_values', eta_x+coeff*sign(unif).*abs(unif).^1.3])';
    else
        if elementBasedSamples
            Xi_values = insertUniform(nurbs.knots{1}, resolution(1));
            Eta_values = insertUniform(nurbs.knots{2}, resolution(2));
        else
            Xi_values = insertUniform3(nurbs.knots{1}, resolution(1), nurbs.degree(1));
            Eta_values = insertUniform3(nurbs.knots{2}, resolution(2), nurbs.degree(2));
        end
    end
%     %% Plot surface where zeta = 1
    X = zeros(length(Xi_values), length(Eta_values));
    Y = zeros(length(Xi_values), length(Eta_values));
    Z = zeros(length(Xi_values), length(Eta_values));
    C = zeros(length(Xi_values), length(Eta_values));
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
    if isnan(d)
%         for i = 1:length(Xi_values)
        parfor i = 1:length(Xi_values)
            X_temp = zeros(1,length(Eta_values));
            Y_temp = zeros(1,length(Eta_values));
            Z_temp = zeros(1,length(Eta_values));
            C_temp = zeros(1,length(Eta_values));
            xi = Xi_values(i);
            for j = 1:length(Eta_values)
                eta = Eta_values(j);
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
                    normal = cross(dydxi,dydeta);
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
                        C_temp(j) = abs(colorFun(y));
                    end
%                     C_temp(j) = colorFun(y.',ny);
%                     y = y + ny*real(C_temp(j))*2;
    %                 C_temp(j)
    %                 C_temp(j) = colorFun(y,xi,eta);
    %                 C_temp(j) = log10(abs(norm(y)-1));
    %                 C_temp(j) = abs(norm(y)-1);
                else
                    y = evaluateNURBS(nurbs, [xi eta]);
                end
                X_temp(j) = y(1);
                Y_temp(j) = y(2);
                Z_temp(j) = y(3);
            end
            X(i,:) = X_temp;
            Y(i,:) = Y_temp;
            Z(i,:) = Z_temp;
            C(i,:) = C_temp;
        end
        if plotSolution
    %         analytic_v = varCol.analytic([reshape(X,length(Xi_values)*length(Eta_values),1), ...
    %                                reshape(Y,length(Xi_values)*length(Eta_values),1), ...
    %                                reshape(Z,length(Xi_values)*length(Eta_values),1)]);
    %         analytic_v = reshape(analytic_v,length(Xi_values),length(Eta_values));
    %         if strcmp(varCol.method,'BEM')
    %             totField = C;
    %             if ~strcmp(varCol.model,'PS')
    %                 P_inc = varCol.P_inc([reshape(X,length(Xi_values)*length(Eta_values),1), ...
    %                                reshape(Y,length(Xi_values)*length(Eta_values),1), ...
    %                                reshape(Z,length(Xi_values)*length(Eta_values),1)]);
    %                 P_inc = reshape(P_inc,length(Xi_values),length(Eta_values));
    %                 C = C - P_inc;
    %             end
    %         else        
    %             if ~strcmp(varCol.model,'PS')
    %                 totField = C + data.P_inc;
    %             end
    %         end
            if false
                noPts = length(Xi_values)*length(Eta_values);
                C = reshape(colorFun([reshape(X,noPts,1),reshape(Y,noPts,1),reshape(Z,noPts,1)],...
                                     [reshape(X,noPts,1),reshape(Y,noPts,1),reshape(Z,noPts,1)].'), ...
                                      length(Xi_values), length(Eta_values));
    %             C = log10(abs(C-analytic_v)/max(max(abs(C))));
            end
            if plotElementEdges
                surf(X,Y,Z, C, 'EdgeColor','none','LineStyle','none')
            else
                surf(X,Y,Z, C)
            end
            maxC = max(max(C));
    %         caxis([max(min(min(C)),-17),max(max(C))])
            colorbar
    %         lighting gouraud
    %         set(gca, 'CScale', 'log')
        else
            if plotElementEdges
                if alphaValue < 1
                    surf(X,Y,Z, 'FaceColor', color,'EdgeColor','none','LineStyle','none','FaceAlpha',alphaValue, 'FaceLighting', 'none')
                else
                    surf(X,Y,Z, 'FaceColor', color,'EdgeColor','none','LineStyle','none','FaceAlpha',alphaValue)
                end
            else
                surf(X,Y,Z, 'FaceColor', color)
            end
        end
        hold on
    %     Xi_values2 = insertUniform3(nurbs.knots{1}, 2, nurbs.degree(1));
    %     Eta_values2 = insertUniform3(nurbs.knots{2}, 2, nurbs.degree(2));
    %     X = zeros(length(Xi_values2), length(Eta_values2));
    %     Y = zeros(length(Xi_values2), length(Eta_values2));
    %     Z = zeros(length(Xi_values2), length(Eta_values2));
    %     Vx = zeros(length(Xi_values2), length(Eta_values2));
    %     Vy = zeros(length(Xi_values2), length(Eta_values2));
    %     Vz = zeros(length(Xi_values2), length(Eta_values2));
    % %     for i = 1:length(Xi_values)
    %     parfor i = 1:length(Xi_values2)
    %         X_temp = zeros(1,length(Eta_values2));
    %         Y_temp = zeros(1,length(Eta_values2));
    %         Z_temp = zeros(1,length(Eta_values2));
    %         V_temp = zeros(3,length(Eta_values2));
    %         xi = Xi_values2(i);
    %         for j = 1:length(Eta_values2)
    %             eta = Eta_values2(j);
    %             [y,dydxi,dydeta] = evaluateNURBS_2ndDeriv(nurbs, [xi eta]);
    %             normal = cross(dydxi,dydeta);
    %             V_temp(:,j) = normal/norm(normal)/2;
    %             X_temp(j) = y(1);
    %             Y_temp(j) = y(2);
    %             Z_temp(j) = y(3);
    %         end
    %         X(i,:) = X_temp;
    %         Y(i,:) = Y_temp;
    %         Z(i,:) = Z_temp;
    %         Vx(i,:) = V_temp(1,:);
    %         Vy(i,:) = V_temp(2,:);
    %         Vz(i,:) = V_temp(3,:);
    %     end
    %     nn = numel(Xi_values2);
    %     mm = numel(Eta_values2);
    %     quiver3(reshape(X,mm*nn,1),reshape(Y,mm*nn,1),reshape(Z,mm*nn,1),reshape(Vx,mm*nn,1),reshape(Vy,mm*nn,1),reshape(Vz,mm*nn,1),0,'color','red')

    else
        nurbspatches{1} = nurbs;
        varCol.dimension = 1;
        varCol = convertNURBS(nurbspatches, varCol);  
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
        hold on
%         parfor e = 1:noElems
        for e = 1:noElems
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


            noXiValues  = max([round(dxi/d)+1, 2]);
            noEtaValues = max([round(deta/d)+1, 2]);
            if noXiValues*noEtaValues < resolution(1)*resolution(2)
                noXiValues = round(resolution(1)*sqrt(dxi/deta));
                noEtaValues = round(resolution(2)*sqrt(deta/dxi));
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
                    surf(X{e},Y{e},Z{e}, C{e}, 'EdgeColor','none','LineStyle','none')
                else
                    surf(X{e},Y{e},Z{e}, C{e})
                end
            else
                if plotElementEdges
                    if alphaValue < 1
                        surf(X{e},Y{e},Z{e}, 'FaceColor', color,'EdgeColor','none','LineStyle','none','FaceAlpha',alphaValue, 'FaceLighting', 'none')
                    else
                        surf(X{e},Y{e},Z{e}, 'FaceColor', color,'EdgeColor','none','LineStyle','none','FaceAlpha',alphaValue)
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
    if plotElementEdges        
        x = zeros(length(Xi_values),1);
        y = zeros(length(Xi_values),1);
        z = zeros(length(Xi_values),1);
        
        for j = 1:length(unique_etaKnots)
            eta = unique_etaKnots(j);
            parfor i = 1:length(Xi_values)
                xi = Xi_values(i);
                v = evaluateNURBS(nurbs, [xi eta]);
                x(i) = v(1);
                y(i) = v(2);
                z(i) = v(3);
            end
            plot3(x,y,z,'color',lineColor,'LineWidth',lineWidth)
        end
        
        x = zeros(length(Eta_values),1);
        y = zeros(length(Eta_values),1);
        z = zeros(length(Eta_values),1);
        for i = 1:length(unique_xiKnots)
            xi = unique_xiKnots(i);
            parfor j = 1:length(Eta_values)
                eta = Eta_values(j);
                v = evaluateNURBS(nurbs, [xi eta]);
                x(j) = v(1);
                y(j) = v(2);
                z(j) = v(3);
            end
            plot3(x,y,z,'color',lineColor,'LineWidth',lineWidth)
        end
    end
% 
%     daspect([1 1 1]);
%     view(3)
%     camlight(110,30);
%     lighting phong;
%     set(gca, 'Color', 'none');
%     camlight('left')
    h = gcf;
    
elseif strcmp(nurbs.type, '2Dsurface')
    
    unique_xiKnots = unique(nurbs.knots{1});
    unique_etaKnots = unique(nurbs.knots{2});
    
    % To reconstruct any interpolatory points, any p repeated knot should be
    % included in the abscissa values to plot
    Xi_values = insertUniform3(nurbs.knots{1}, resolution(1), nurbs.degree(1));
    Eta_values = insertUniform3(nurbs.knots{2}, resolution(2), nurbs.degree(2));

    X = zeros(2*length(Xi_values)+2*length(Eta_values)-3,1);
    Y = zeros(2*length(Xi_values)+2*length(Eta_values)-3,1);
    counter = 1;
    % along eta = 0
    for i = 1:length(Xi_values)
        xi = Xi_values(i);
        eta = 0;
        v = evaluateNURBS(nurbs, [xi eta]);
        X(counter) = v(1);
        Y(counter) = v(2);
        counter = counter + 1;
    end
    % along xi = 1
    for i = 2:length(Eta_values)
        xi = 1;
        eta = Eta_values(i);
        v = evaluateNURBS(nurbs, [xi eta]);
        X(counter) = v(1);
        Y(counter) = v(2);
        counter = counter + 1;
    end
    % reverse order of arrays
    Xi_values = Xi_values(end:-1:1);
    Eta_values = Eta_values(end:-1:1);
    
    % along eta = 1
    for i = 2:length(Xi_values)
        xi = Xi_values(i);
        eta = 1;
        v = evaluateNURBS(nurbs, [xi eta]);
        X(counter) = v(1);
        Y(counter) = v(2);
        counter = counter + 1;
    end
    % along xi = 0
    for i = 2:length(Eta_values)
        xi = 0;
        eta = Eta_values(i);
        v = evaluateNURBS(nurbs, [xi eta]);
        X(counter) = v(1);
        Y(counter) = v(2);
        counter = counter + 1;
    end    
    % reverse back order of arrays
    Xi_values = Xi_values(end:-1:1);
    Eta_values = Eta_values(end:-1:1);
    if plotElementEdges
        fill(X,Y, color,'EdgeColor','none','LineStyle','none')
    else
        fill(X,Y, color)
    end
    hold on

    if plotElementEdges        
        % Plot element mesh for constant xi
        x = zeros(length(Eta_values),1);
        y = zeros(length(Eta_values),1);
        
        for i = 1:length(unique_xiKnots)
            xi = unique_xiKnots(i);
            for j = 1:length(Eta_values)
                eta = Eta_values(j);
                v = evaluateNURBS(nurbs, [xi eta]);
                x(j) = v(1);
                y(j) = v(2);
            end
            plot(x,y,'color',lineColor)
        end

        % Plot element mesh for constant eta
        x = zeros(length(Xi_values),1);
        y = zeros(length(Xi_values),1);
        
        for j = 1:length(unique_etaKnots)
            eta = unique_etaKnots(j);
            for i = 1:length(Xi_values)
                xi = Xi_values(i);
                v = evaluateNURBS(nurbs, [xi eta]);
                x(i) = v(1);
                y(i) = v(2);
            end
            plot(x,y,'color',lineColor)
        end
        
    end
    
elseif strcmp(nurbs.type, '3Dcurve')
    Xi_values = insertUniform(nurbs.knots, resolution);

    C = zeros(length(Xi_values), 3);

    for j = 1:length(Xi_values)
        xi = Xi_values(j);
        C(j,:) = evaluateNURBS(nurbs, xi);
    end

    h = plot3(C(:,1), C(:,2), C(:,3), 'color', color); 
elseif strcmp(nurbs.type, '2Dcurve')
    Xi_values = insertUniform(nurbs.knots, resolution);

    C = zeros(length(Xi_values), 2);

    for j = 1:length(Xi_values)
        xi = Xi_values(j);
        C(j,:) = evaluateNURBS(nurbs, xi);
    end

    h = plot(C(:,1), C(:,2), 'color', color);  
elseif strcmp(nurbs.type, '1Dnurbs')
    Xi_values = insertUniform(nurbs.knots, resolution);

    C = zeros(length(Xi_values), 1);

    for j = 1:length(Xi_values)
        xi = Xi_values(j);
        C(j) = evaluateNURBS(nurbs, xi);
    end

    h = plot(Xi_values, C, 'color', color);    
end






