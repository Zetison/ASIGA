function plotLRinParameterSpace(LRobj, resolution, basisNrXi, basisNrEta)
if nargin < 3
    plotElementEdges = true;
end
if nargin < 4
    % Plots the surface of a NURBS volume
    % colors = colormap('summer');
    mycolor = [33 76 161]/255;
end
color = 1.5*[44 77 32]/255; % colors(end,:);
if nargin < 5
    alphaValue = 1;
end
% color = [33 76 161]/255; %214CA1
lineColor = 'black';
if strcmp(LRobj.type, '3Dvolume')
    
elseif strcmp(LRobj.type, '3Dsurface')
   
elseif strcmp(LRobj.type, '2Dsurface')
    Xi = LRobj.knots{1};
    Eta = LRobj.knots{2};
    n_xi = LRobj.number(1);
    n_eta = LRobj.number(2);
    p_xi = LRobj.degree(1);
    p_eta = LRobj.degree(2);
    
    unique_xiKnots = unique(Xi);
    unique_etaKnots = unique(Eta);
    
    % To reconstruct any interpolatory points, any p repeated knot should be
    % included in the abscissa values to plot
    Xi_values = insertUniform3(LRobj.knots{1}, resolution(1), LRobj.degree(1));
    Eta_values = insertUniform3(LRobj.knots{2}, resolution(2), LRobj.degree(2));

    X = zeros(length(Xi_values),length(Eta_values));
    Y = X;
    Z = X;
    for i = 1:length(Xi_values)
        xi = Xi_values(i);
        for j = 1:length(Eta_values)
            eta = Eta_values(j);
            X(i,j) = xi;
            Y(i,j) = eta;
            
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

    M_xi = LRobj.M_xi;
    M_eta = LRobj.M_eta;
    
    % Plot element mesh for constant xi
    x = zeros(length(Eta_values),1);
    y = zeros(length(Eta_values),1);
    M_xi(counter).eta = uniqueEta(i);
    M_xi(counter).xiRange = [Xi(1) Xi(end)];
            
    for i = 1:length(M_xi)
        xi = M_xi(i).xi;
        z = zeros(length(Eta_values),1);
        Eta_values = linspace(M_xi(i).etaRange(1),M_xi(i).etaRange(2),resolution(2));
        for j = 1:length(Eta_values)
            eta = Eta_values(j);
            x(j) = xi;
            y(j) = eta;
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

    for i = 1:length(M_eta)
        eta = M_eta(i).eta;
        z = zeros(length(Xi_values),1);
        Xi_values = linspace(M_eta(i).xiRange(1),M_eta(i).xiRange(2),resolution(1));
        for j = 1:length(Xi_values)
            xi = Xi_values(j);
            x(j) = xi;
            y(j) = eta;
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
        
elseif strcmp(LRobj.type, '2Dcurve')  
end






