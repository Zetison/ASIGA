function h = plotCuttingPlaneConstEta(nurbs, resolution, plane, eta, lineColor)
if nargin < 5
    lineColor = 'black';
end

% color = [33 76 161]/255; %214CA1
if strcmp(nurbs.type, '3Dvolume')

    unique_xiKnots = unique(nurbs.knots{1});
    unique_zetaKnots = unique(nurbs.knots{3});
    
    % To reconstruct any interpolatory points, any p repeated knot should be
    % included in the abscissa values to plot
    Xi_values = insertUniform3(nurbs.knots{1}, resolution(1), 1);
    Zeta_values = insertUniform3(nurbs.knots{3}, resolution(2), 1);

    x = zeros(length(Xi_values),1);
    y = zeros(length(Xi_values),1);
    for k = 1:length(unique_zetaKnots)
        zeta = unique_zetaKnots(k);
        parfor j = 1:length(Xi_values)
            xi = Xi_values(j);
            v = evaluateNURBS(nurbs, [xi eta zeta]);
            switch plane
                case 'xz'
                    x(j) = v(1);
                    y(j) = v(3);
                case 'yz'
                    x(j) = v(2);
                    y(j) = v(3);
            end
        end
        plot(x,y,'color',lineColor)
        hold on
    end

    x = zeros(length(Zeta_values),1);
    y = zeros(length(Zeta_values),1);
    for j = 1:length(unique_xiKnots)
        xi = unique_xiKnots(j);
        parfor k = 1:length(Zeta_values)
            zeta = Zeta_values(k);
            v = evaluateNURBS(nurbs, [xi eta zeta]);
            switch plane
                case 'xz'
                    x(k) = v(1);
                    y(k) = v(3);
                case 'yz'
                    x(k) = v(2);
                    y(k) = v(3);
            end
        end
        plot(x,y,'color',lineColor)
    end
% 
%     daspect([1 1 1]);
%     view(3)
%     camlight(110,30);
%     lighting phong;
%     set(gca, 'Color', 'none');
%     camlight('left')
    h = gcf;
else
    error('Not implemented!')
end






