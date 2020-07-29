function h = plotCuttingPlaneConstXi(nurbs, resolution, plane, xi, lineColor)
error('Depricated. Use plotNURBS() with cutting planes at the end of patches')
if nargin < 5
    lineColor = 'black';
end

% color = [33 76 161]/255; %214CA1
if strcmp(nurbs.type, '3Dvolume')

    unique_etaKnots = unique(nurbs.knots{2});
    unique_zetaKnots = unique(nurbs.knots{3});
    
    % To reconstruct any interpolatory points, any p repeated knot should be
    % included in the abscissa values to plot
    Eta_values = insertUniform3(nurbs.knots{2}, resolution(1), 1);
    Zeta_values = insertUniform3(nurbs.knots{3}, resolution(2), 1);

    x = zeros(length(Eta_values),length(unique_zetaKnots));
    y = zeros(length(Eta_values),length(unique_zetaKnots));
    parfor k = 1:length(unique_zetaKnots)
        zeta = unique_zetaKnots(k);
        x_temp = zeros(length(Eta_values),1);
        y_temp = zeros(length(Eta_values),1);
        for j = 1:length(Eta_values)
            eta = Eta_values(j);
            v = evaluateNURBS(nurbs, [xi eta zeta]);
            switch plane
                case 'xy'
                    x_temp(j) = v(1);
                    y_temp(j) = v(2);
                case 'xz'
                    x_temp(j) = v(1);
                    y_temp(j) = v(3);
            end
        end
        x(:,k) = x_temp;
        y(:,k) = y_temp;
    end
    for k = 1:length(unique_zetaKnots)
        plot(x(:,k),y(:,k),'color',lineColor)
        hold on
    end
        

    x = zeros(length(Zeta_values),length(unique_etaKnots));
    y = zeros(length(Zeta_values),length(unique_etaKnots));
    parfor j = 1:length(unique_etaKnots)
        eta = unique_etaKnots(j);
        x_temp = zeros(length(Zeta_values),1);
        y_temp = zeros(length(Zeta_values),1);
        for k = 1:length(Zeta_values)
            zeta = Zeta_values(k);
            v = evaluateNURBS(nurbs, [xi eta zeta]);
            switch plane
                case 'xy'
                    x_temp(k) = v(1);
                    y_temp(k) = v(2);
                case 'xz'
                    x_temp(k) = v(1);
                    y_temp(k) = v(3);
            end
        end
        x(:,j) = x_temp;
        y(:,j) = y_temp;
    end
    for k = 1:length(unique_etaKnots)
        plot(x(:,k),y(:,k),'color',lineColor)
        hold on
    end
    
    h = gcf;
else
    error('Not implemented!')
end






