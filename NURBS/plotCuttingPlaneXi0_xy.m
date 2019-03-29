function h = plotCuttingPlaneXi0_xy(nurbs, resolution)

% color = [33 76 161]/255; %214CA1
lineColor = 'black';
if strcmp(nurbs.type, '3Dvolume')

    unique_etaKnots = unique(nurbs.knots{2});
    unique_zetaKnots = unique(nurbs.knots{3});
    xi = 0;
    % To reconstruct any interpolatory points, any p repeated knot should be
    % included in the abscissa values to plot
    Eta_values = insertUniform3(nurbs.knots{2}, resolution(1), 1);
    Zeta_values = insertUniform3(nurbs.knots{3}, resolution(2), 1);

    x = zeros(length(Eta_values),1);
    y = zeros(length(Eta_values),1);
    for k = 1:length(unique_zetaKnots)
        zeta = unique_zetaKnots(k);
        for j = 1:length(Eta_values)
            eta = Eta_values(j);
            v = evaluateNURBS(nurbs, [xi eta zeta]);
            x(j) = v(1);
            y(j) = v(2);
        end
        plot(x,y,'color',lineColor)
        hold on
    end

    x = zeros(length(Zeta_values),1);
    y = zeros(length(Zeta_values),1);
    for j = 1:length(unique_etaKnots)
        eta = unique_etaKnots(j);
        for k = 1:length(Zeta_values)
            zeta = Zeta_values(k);
            v = evaluateNURBS(nurbs, [xi eta zeta]);
            x(k) = v(1);
            y(k) = v(2);
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






