function plotNURBSparametricSpace(nurbs, resolution, basisNrXi, basisNrEta)
error('Depricated...')
if nargin < 3
    plotElementEdges = true;
end
if nargin < 4
    % Plots the surface of a NURBS volume
    % colors = colormap('summer');
    mycolor = [33 76 161]/255;
end
color = getColor(1); % colors(end,:);
if nargin < 5
    alphaValue = 1;
end
% color = [33 76 161]/255; %214CA1
lineColor = 'black';
if strcmp(nurbs.type, '3Dvolume')
    
elseif strcmp(nurbs.type, '3Dsurface')
   
elseif strcmp(nurbs.type, '2Dsurface')
    Xi = nurbs.knots{1};
    Eta = nurbs.knots{2};
    n_xi = nurbs.number(1);
    n_eta = nurbs.number(2);
    p_xi = nurbs.degree(1);
    p_eta = nurbs.degree(2);
    
    unique_xiKnots = unique(Xi);
    unique_etaKnots = unique(Eta);
    
    % To reconstruct any interpolatory points, any p repeated knot should be
    % included in the abscissa values to plot
    Xi_values = insertUniform3(nurbs.knots{1}, resolution(1), nurbs.degree(1));
    Eta_values = insertUniform3(nurbs.knots{2}, resolution(2), nurbs.degree(2));

    X = zeros(length(Xi_values),length(Eta_values));
    Y = X;
    Z = X;
    for i = 1:length(Xi_values)
        xi = Xi_values(i);
        for j = 1:length(Eta_values)
            eta = Eta_values(j);
            v = evaluateNURBS(nurbs, [xi eta]);
            X(i,j) = v(1);
            Y(i,j) = v(2);
            
            if ((Xi(basisNrXi) <= xi && xi < Xi(basisNrXi+p_xi+1)) || (xi == Xi(end) && basisNrXi == n_xi)) && ...
               ((Eta(basisNrXi) <= eta && eta < Eta(basisNrEta+p_eta+1)) || (eta == Eta(end) && basisNrEta == n_eta))
                i1 = findKnotSpan(n_xi, p_xi, xi, Xi);
                i2 = findKnotSpan(n_eta, p_eta, eta, Eta);
                N = Bspline_basis(i1, xi, p_xi, Xi, 0);
                M = Bspline_basis(i2, eta, p_eta, Eta, 0);

                Z(i,j) = N(p_xi+1+basisNrXi-i1)*M(p_eta+1+basisNrEta-i2);
            end
        end
    end
    % reverse back order of arrays
    Xi_values = Xi_values(end:-1:1);
    Eta_values = Eta_values(end:-1:1);
    surf(X,Y, Z,'EdgeColor','none','LineStyle','none')
    customColorMap  = getCustomColorMap(); 
    customColorMap = customColorMap(((end+1)/2):end,:);
    colormap(customColorMap)
    camlight
    axis off
    hold on

    % Plot element mesh for constant xi
    x = zeros(length(Eta_values),1);
    y = zeros(length(Eta_values),1);

    for i = 1:length(unique_xiKnots)
        xi = unique_xiKnots(i);
        z = zeros(length(Eta_values),1);
        for j = 1:length(Eta_values)
            eta = Eta_values(j);
            v = evaluateNURBS(nurbs, [xi eta]);
            x(j) = v(1);
            y(j) = v(2);
            if ((Xi(basisNrXi) <= xi && xi < Xi(basisNrXi+p_xi+1)) || (xi == Xi(end) && basisNrXi == n_xi)) && ...
               ((Eta(basisNrXi) <= eta && eta < Eta(basisNrEta+p_eta+1)) || (eta == Eta(end) && basisNrEta == n_eta))
                i1 = findKnotSpan(n_xi, p_xi, xi, Xi);
                i2 = findKnotSpan(n_eta, p_eta, eta, Eta);
                N = Bspline_basis(i1, xi, p_xi, Xi, 0);
                M = Bspline_basis(i2, eta, p_eta, Eta, 0);

                z(j) = N(p_xi+1+basisNrXi-i1)*M(p_eta+1+basisNrEta-i2);
            end
        end
        plot3(x,y,z,'color',lineColor)
    end

    % Plot element mesh for constant eta
    x = zeros(length(Xi_values),1);
    y = zeros(length(Xi_values),1);

    for j = 1:length(unique_etaKnots)
        eta = unique_etaKnots(j);
        z = zeros(length(Xi_values),1);
        for i = 1:length(Xi_values)
            xi = Xi_values(i);
            v = evaluateNURBS(nurbs, [xi eta]);
            x(i) = v(1);
            y(i) = v(2);
            if ((Xi(basisNrXi) <= xi && xi < Xi(basisNrXi+p_xi+1)) || (xi == Xi(end) && basisNrXi == n_xi)) && ...
               ((Eta(basisNrXi) <= eta && eta < Eta(basisNrEta+p_eta+1)) || (eta == Eta(end) && basisNrEta == n_eta))
                i1 = findKnotSpan(n_xi, p_xi, xi, Xi);
                i2 = findKnotSpan(n_eta, p_eta, eta, Eta);
                N = Bspline_basis(i1, xi, p_xi, Xi, 0);
                M = Bspline_basis(i2, eta, p_eta, Eta, 0);

                z(i) = N(p_xi+1+basisNrXi-i1)*M(p_eta+1+basisNrEta-i2);
            end
        end
        plot3(x,y,z,'color',lineColor)
    end
        
elseif strcmp(nurbs.type, '2Dcurve')  
end






