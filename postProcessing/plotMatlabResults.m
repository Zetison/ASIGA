function [X, Y, W] = plotResults(nurbs, results, scalarfield, controlPts, weights, resolution, plotEdge, mirrorMatrix, ...
    plotElementEdges, color)

signFactor = 1; % set to -1 to reverse axis direction
customColorMap  = getCustomColorMap(); 

if nargin < 7
    plotEdge = [1 1;    % At xi = 0 and xi = 1
                1 1;    % At eta = 0 and eta = 1
                1 1];   % At zeta = 0 and zeta = 1
end
if nargin < 8
    mirrorMatrix = [];
end
if nargin < 9
    plotElementEdges = true;
end
if nargin < 10
    % Plots the surface of a NURBS volume
    % colors = colormap('summer');
    color = 1.5*[44 77 32]/255; % colors(end,:);
    mycolor = [33 76 161]/255;
end
% color = [33 76 161]/255; %214CA1
lineColor = 'black';


if strcmp(nurbs.type, '3Dvolume')
    maxW = -inf;
    minW = inf;
    p = nurbs.degree(1);
    q = nurbs.degree(2);
    r = nurbs.degree(3);
    
    n = nurbs.number(1);
    m = nurbs.number(2);
    l = nurbs.number(3);
    
    Xi = nurbs.knots{1};
    Eta = nurbs.knots{2};
    Zeta = nurbs.knots{3};
    
    Ux = results{1};
    Uy = results{2};
    Uz = results{3};
    
    noElemsXi       = length(unique(Xi))-1; % # of elements xi dir.
    noElemsEta       = length(unique(Eta))-1; % # of elements xi dir.
    
    unique_xiKnots = unique(nurbs.knots{1});
    unique_etaKnots = unique(nurbs.knots{2});
    unique_zetaKnots = unique(nurbs.knots{3});
    
    % To reconstruct any interpolatory points, any p repeated knot should be
    % included in the abscissa values to plot
    Xi_values = insertUniform3(nurbs.knots{1}, resolution(1), nurbs.degree(1));
    Eta_values = insertUniform3(nurbs.knots{2}, resolution(2), nurbs.degree(2));
    Zeta_values = insertUniform3(nurbs.knots{3}, resolution(3), nurbs.degree(3));
    
    if plotEdge(1,1)
        %% Plot surface where xi = 0
        xi0Nodes = zeros(1,m*l);
        counter = 1;
        for k = 1:l
            for j = 1:m
                xi0Nodes(counter) = (m*n)*(k-1) + n*(j-1) + 1;
                counter = counter + 1;
            end
        end

        [EtaZetaMesh, ~, ~, elRangeEta, elRangeZeta] = generateIGA2DMesh(Eta, Zeta, q, r, m, l);
 
        XYZ = zeros(length(Eta_values), length(Zeta_values),3);
        W = zeros(length(Eta_values), length(Zeta_values));
        for j = 1:length(Eta_values)
            for k = 1:length(Zeta_values)
                eta = Eta_values(j);
                zeta = Zeta_values(k);
                [R_fun, ~, ~] = NURBS2DBasis(eta, zeta, q, r, Eta, Zeta, weights(xi0Nodes));


                idx_eta = findElementInMesh(eta, elRangeEta);
                idx_zeta = findElementInMesh(zeta, elRangeZeta);

                e = noElemsEta*(idx_zeta-1) + idx_eta;

                sctrEtaZeta = xi0Nodes(EtaZetaMesh(e,:));          %  element scatter vector
                v = R_fun*controlPts(sctrEtaZeta,:);
                u = R_fun*[Ux(sctrEtaZeta), Uy(sctrEtaZeta), Uz(sctrEtaZeta)];
                
                switch scalarfield
                    case 1
                        W(j,k) = signFactor*u(1);
                    case 2
                        W(j,k) = signFactor*u(2);
                    case 3
                        W(j,k) = signFactor*u(3);
                    case 4
                        W(j,k) = norm(u);
                    case 5
                        x = v(1);
                        y = v(2);
                        u_x = signFactor*u(1);
                        u_y = signFactor*u(2);
                        theta = atan2(y,x);
                        W(j,k) = u_y*sin(theta) + u_x*cos(theta);
                    case 6
                        x = v(1);
                        y = v(2);
                        z = v(3);
                        radius = sqrt(x^2+y^2+z^2);
                        phi = atan2(y,x);
                        theta = acos(z/radius);

                        u_x = signFactor*u(1);
                        u_y = signFactor*u(2);
                        u_z = signFactor*u(3);

                        W(j,k) = sin(theta)*cos(phi)*u_x + sin(theta)*sin(phi)*u_y + cos(theta)*u_z;
                end
                XYZ(j,k,:) = u + v;
            end
        end
        if minW > min(min(W))
            minW = min(min(W));
        end
        if maxW < max(max(W))
            maxW = max(max(W));
        end
        plotXYZwithMirroring
        
        if plotElementEdges
            % Plot element mesh where xi = 0
            xyz = zeros(3, length(Zeta_values));
            for j = 1:length(unique_etaKnots)
                eta = unique_etaKnots(j);
                for k = 1:length(Zeta_values)
                    zeta = Zeta_values(k);
                    [R_fun, ~, ~] = NURBS2DBasis(eta, zeta, q, r, Eta, Zeta, weights(xi0Nodes));


                    idx_eta = findElementInMesh(eta, elRangeEta);
                    idx_zeta = findElementInMesh(zeta, elRangeZeta);

                    e = noElemsEta*(idx_zeta-1) + idx_eta;

                    sctrEtaZeta = xi0Nodes(EtaZetaMesh(e,:));          %  element scatter vector
                    v = R_fun*controlPts(sctrEtaZeta,:);
                    u = R_fun*[Ux(sctrEtaZeta), Uy(sctrEtaZeta), Uz(sctrEtaZeta)];
                    xyz(:,k) = u + v;
                end
                plot_xyz_withMirroring
            end

            xyz = zeros(3, length(Eta_values));
            for k = 1:length(unique_zetaKnots)
                zeta = unique_zetaKnots(k);
                for j = 1:length(Eta_values)
                    [R_fun, ~, ~] = NURBS2DBasis(eta, zeta, q, r, Eta, Zeta, weights(xi0Nodes));


                    idx_eta = findElementInMesh(eta, elRangeEta);
                    idx_zeta = findElementInMesh(zeta, elRangeZeta);

                    e = noElemsEta*(idx_zeta-1) + idx_eta;

                    sctrEtaZeta = xi0Nodes(EtaZetaMesh(e,:));          %  element scatter vector
                    v = R_fun*controlPts(sctrEtaZeta,:);
                    u = R_fun*[Ux(sctrEtaZeta), Uy(sctrEtaZeta), Uz(sctrEtaZeta)];
                    xyz(:,j) = u + v;
                end
                plot_xyz_withMirroring
            end
        end
    end

    if plotEdge(1,2)
        %% Plot surface where xi = 1
        xi1Nodes = zeros(1,m*l);
        counter = 1;
        for k = 1:l
            for j = 1:m
                xi1Nodes(counter) = (m*n)*(k-1) + n*(j-1) + n;
                counter = counter + 1;
            end
        end
        [EtaZetaMesh, ~, ~, elRangeEta, elRangeZeta]  = generateIGA2DMesh(Eta, Zeta, q, r, m, l);
     
        XYZ = zeros(length(Eta_values), length(Zeta_values),3);
        W = zeros(length(Eta_values), length(Zeta_values));
        for j = 1:length(Eta_values)
            for k = 1:length(Zeta_values)
                eta = Eta_values(j);
                zeta = Zeta_values(k);
                [R_fun, ~, ~] = NURBS2DBasis(eta, zeta, q, r, Eta, Zeta, weights(xi1Nodes));


                idx_eta = findElementInMesh(eta, elRangeEta);
                idx_zeta = findElementInMesh(zeta, elRangeZeta);

                e = noElemsEta*(idx_zeta-1) + idx_eta;

                sctrEtaZeta = xi1Nodes(EtaZetaMesh(e,:));          %  element scatter vector
                v = R_fun*controlPts(sctrEtaZeta,:);
                u = R_fun*[Ux(sctrEtaZeta), Uy(sctrEtaZeta), Uz(sctrEtaZeta)];
                switch scalarfield
                    case 1
                        W(j,k) = signFactor*u(1);
                    case 2
                        W(j,k) = signFactor*u(2);
                    case 3
                        W(j,k) = signFactor*u(3);
                    case 4
                        W(j,k) = norm(u);
                    case 5
                        x = v(1);
                        y = v(2);
                        u_x = signFactor*u(1);
                        u_y = signFactor*u(2);
                        theta = atan2(y,x);
                        W(j,k) = u_y*sin(theta) + u_x*cos(theta);
                    case 6
                        x = v(1);
                        y = v(2);
                        z = v(3);
                        radius = sqrt(x^2+y^2+z^2);
                        phi = atan2(y,x);
                        theta = acos(z/radius);

                        u_x = signFactor*u(1);
                        u_y = signFactor*u(2);
                        u_z = signFactor*u(3);

                        W(j,k) = sin(theta)*cos(phi)*u_x + sin(theta)*sin(phi)*u_y + cos(theta)*u_z;
                end
                XYZ(j,k,:) = u + v;
            end
        end
        if minW > min(min(W))
            minW = min(min(W));
        end
        if maxW < max(max(W))
            maxW = max(max(W));
        end
        plotXYZwithMirroring
        
        if plotElementEdges
            % Plot element mesh where xi = 1
            xyz = zeros(3, length(Zeta_values));
            for j = 1:length(unique_etaKnots)
                eta = unique_etaKnots(j);
                for k = 1:length(Zeta_values)
                    zeta = Zeta_values(k);
                    [R_fun, ~, ~] = NURBS2DBasis(eta, zeta, q, r, Eta, Zeta, weights(xi1Nodes));


                    idx_eta = findElementInMesh(eta, elRangeEta);
                    idx_zeta = findElementInMesh(zeta, elRangeZeta);

                    e = noElemsEta*(idx_zeta-1) + idx_eta;

                    sctrEtaZeta = xi1Nodes(EtaZetaMesh(e,:));          %  element scatter vector
                    v = R_fun*controlPts(sctrEtaZeta,:);
                    u = R_fun*[Ux(sctrEtaZeta), Uy(sctrEtaZeta), Uz(sctrEtaZeta)];
                    xyz(:,k) = u + v;
                end
                plot_xyz_withMirroring
            end

            xyz = zeros(3, length(Eta_values));
            for k = 1:length(unique_zetaKnots)
                zeta = unique_zetaKnots(k);
                for j = 1:length(Eta_values)
                    eta = Eta_values(j);
                    [R_fun, ~, ~] = NURBS2DBasis(eta, zeta, q, r, Eta, Zeta, weights(xi1Nodes));


                    idx_eta = findElementInMesh(eta, elRangeEta);
                    idx_zeta = findElementInMesh(zeta, elRangeZeta);

                    e = noElemsEta*(idx_zeta-1) + idx_eta;

                    sctrEtaZeta = xi1Nodes(EtaZetaMesh(e,:));          %  element scatter vector
                    v = R_fun*controlPts(sctrEtaZeta,:);
                    u = R_fun*[Ux(sctrEtaZeta), Uy(sctrEtaZeta), Uz(sctrEtaZeta)];
                    xyz(:,j) = u + v;
                end
                plot_xyz_withMirroring
            end
        end
    end
    
    if plotEdge(2,1)
        %% Plot surface where eta = 0
        eta0Nodes = zeros(1,n*l);
        counter = 1;
        for k = 1:l
            for i = 1:n
                eta0Nodes(counter) = (m*n)*(k-1) + i;
                counter = counter + 1;
            end
        end

        [XiZetaMesh, ~, ~, elRangeXi, elRangeZeta]  = generateIGA2DMesh(Xi, Zeta, p, r, n, l);
     
        XYZ = zeros(length(Xi_values), length(Zeta_values),3);
        W = zeros(length(Xi_values), length(Zeta_values));
        for i = 1:length(Xi_values)
            for k = 1:length(Zeta_values)
                xi = Xi_values(i);
                zeta = Zeta_values(k);
                [R_fun, ~, ~] = NURBS2DBasis(xi, zeta, p, r, Xi, Zeta, weights(eta0Nodes));


                idx_xi = findElementInMesh(xi, elRangeXi);
                idx_zeta = findElementInMesh(zeta, elRangeZeta);

                e = noElemsXi*(idx_zeta-1) + idx_xi;

                sctrXiZeta = eta0Nodes(XiZetaMesh(e,:));          %  element scatter vector
                v = R_fun*controlPts(sctrXiZeta,:);
                u = R_fun*[Ux(sctrXiZeta), Uy(sctrXiZeta), Uz(sctrXiZeta)];
                switch scalarfield
                    case 1
                        W(i,k) = signFactor*u(1);
                    case 2
                        W(i,k) = signFactor*u(2);
                    case 3
                        W(i,k) = signFactor*u(3);
                    case 4
                        W(i,k) = norm(u);
                    case 5
                        x = v(1);
                        y = v(2);
                        u_x = signFactor*u(1);
                        u_y = signFactor*u(2);
                        theta = atan2(y,x);
                        W(i,k) = u_y*sin(theta) + u_x*cos(theta);
                    case 6
                        x = v(1);
                        y = v(2);
                        z = v(3);
                        radius = sqrt(x^2+y^2+z^2);
                        phi = atan2(y,x);
                        theta = acos(z/radius);

                        u_x = signFactor*u(1);
                        u_y = signFactor*u(2);
                        u_z = signFactor*u(3);

                        W(i,k) = sin(theta)*cos(phi)*u_x + sin(theta)*sin(phi)*u_y + cos(theta)*u_z;
                end
                XYZ(i,k,:) = u + v;
            end
        end
        if minW > min(min(W))
            minW = min(min(W));
        end
        if maxW < max(max(W))
            maxW = max(max(W));
        end
        plotXYZwithMirroring
        
        if plotElementEdges
            % Plot element mesh where eta = 0
            xyz = zeros(3, length(Zeta_values));
            for i = 1:length(unique_xiKnots)
                xi = unique_xiKnots(i);
                for k = 1:length(Zeta_values)
                    zeta = Zeta_values(k);
                    [R_fun, ~, ~] = NURBS2DBasis(xi, zeta, p, r, Xi, Zeta, weights(eta0Nodes));


                    idx_xi = findElementInMesh(xi, elRangeXi);
                    idx_zeta = findElementInMesh(zeta, elRangeZeta);

                    e = noElemsXi*(idx_zeta-1) + idx_xi;

                    sctrXiZeta = eta0Nodes(XiZetaMesh(e,:));          %  element scatter vector
                    v = R_fun*controlPts(sctrXiZeta,:);
                    u = R_fun*[Ux(sctrXiZeta), Uy(sctrXiZeta), Uz(sctrXiZeta)];
                    xyz(:,k) = u + v;
                end
                plot_xyz_withMirroring
            end

            xyz = zeros(3, length(Xi_values));
            for k = 1:length(unique_zetaKnots)
                zeta = unique_zetaKnots(k);
                for i = 1:length(Xi_values)
                    xi = Xi_values(i);
                    [R_fun, ~, ~] = NURBS2DBasis(xi, zeta, p, r, Xi, Zeta, weights(eta0Nodes));


                    idx_xi = findElementInMesh(xi, elRangeXi);
                    idx_zeta = findElementInMesh(zeta, elRangeZeta);

                    e = noElemsXi*(idx_zeta-1) + idx_xi;

                    sctrXiZeta = eta0Nodes(XiZetaMesh(e,:));          %  element scatter vector
                    v = R_fun*controlPts(sctrXiZeta,:);
                    u = R_fun*[Ux(sctrXiZeta), Uy(sctrXiZeta), Uz(sctrXiZeta)];
                    xyz(:,i) = u + v;
                end
                plot_xyz_withMirroring
            end
        end
    end

    if plotEdge(2,2)
        %% Plot surface where eta = 1
        eta1Nodes = zeros(1,n*l);
        counter = 1;
        for k = 1:l
            for i = 1:n
                eta1Nodes(counter) = (m*n)*(k-1) + n*(m-1) + i;
                counter = counter + 1;
            end
        end

        [XiZetaMesh, ~, ~, elRangeXi, elRangeZeta]  = generateIGA2DMesh(Xi, Zeta, p, r, n, l);
     
        XYZ = zeros(length(Xi_values), length(Zeta_values),3);
        W = zeros(length(Xi_values), length(Zeta_values));
        for i = 1:length(Xi_values)
            for k = 1:length(Zeta_values)
                xi = Xi_values(i);
                zeta = Zeta_values(k);
                [R_fun, ~, ~] = NURBS2DBasis(xi, zeta, p, r, Xi, Zeta, weights(eta1Nodes));


                idx_xi = findElementInMesh(xi, elRangeXi);
                idx_zeta = findElementInMesh(zeta, elRangeZeta);

                e = noElemsXi*(idx_zeta-1) + idx_xi;

                sctrXiZeta = eta1Nodes(XiZetaMesh(e,:));          %  element scatter vector
                v = R_fun*controlPts(sctrXiZeta,:);
                u = R_fun*[Ux(sctrXiZeta), Uy(sctrXiZeta), Uz(sctrXiZeta)];
                switch scalarfield
                    case 1
                        W(i,k) = signFactor*u(1);
                    case 2
                        W(i,k) = signFactor*u(2);
                    case 3
                        W(i,k) = signFactor*u(3);
                    case 4
                        W(i,k) = norm(u);
                    case 5
                        x = v(1);
                        y = v(2);
                        u_x = signFactor*u(1);
                        u_y = signFactor*u(2);
                        theta = atan2(y,x);
                        W(i,k) = u_y*sin(theta) + u_x*cos(theta);
                    case 6
                        x = v(1);
                        y = v(2);
                        z = v(3);
                        radius = sqrt(x^2+y^2+z^2);
                        phi = atan2(y,x);
                        theta = acos(z/radius);

                        u_x = signFactor*u(1);
                        u_y = signFactor*u(2);
                        u_z = signFactor*u(3);

                        W(i,k) = sin(theta)*cos(phi)*u_x + sin(theta)*sin(phi)*u_y + cos(theta)*u_z;
                end
                XYZ(i,k,:) = u + v;
            end
        end
        if minW > min(min(W))
            minW = min(min(W));
        end
        if maxW < max(max(W))
            maxW = max(max(W));
        end
        plotXYZwithMirroring
        
        if plotElementEdges
            % Plot element mesh where eta = 1
            xyz = zeros(3, length(Zeta_values));
            for i = 1:length(unique_xiKnots)
                xi = unique_xiKnots(i);
                for k = 1:length(Zeta_values)
                    zeta = Zeta_values(k);
                    [R_fun, ~, ~] = NURBS2DBasis(xi, zeta, p, r, Xi, Zeta, weights(eta1Nodes));

                    idx_xi = findElementInMesh(xi, elRangeXi);
                    idx_zeta = findElementInMesh(zeta, elRangeZeta);

                    e = noElemsXi*(idx_zeta-1) + idx_xi;

                    sctrXiZeta = eta1Nodes(XiZetaMesh(e,:));          %  element scatter vector
                    v = R_fun*controlPts(sctrXiZeta,:);
                    u = R_fun*[Ux(sctrXiZeta), Uy(sctrXiZeta), Uz(sctrXiZeta)];                    
                    
                    xyz(:,k) = u + v;
                end
                plot_xyz_withMirroring
            end

            xyz = zeros(3, length(Xi_values));
            for k = 1:length(unique_zetaKnots)
                zeta = unique_zetaKnots(k);
                for i = 1:length(Xi_values)
                    xi = Xi_values(i);
                    [R_fun, ~, ~] = NURBS2DBasis(xi, zeta, p, r, Xi, Zeta, weights(eta1Nodes));

                    idx_xi = findElementInMesh(xi, elRangeXi);
                    idx_zeta = findElementInMesh(zeta, elRangeZeta);

                    e = noElemsXi*(idx_zeta-1) + idx_xi;

                    sctrXiZeta = eta1Nodes(XiZetaMesh(e,:));          %  element scatter vector
                    v = R_fun*controlPts(sctrXiZeta,:);
                    u = R_fun*[Ux(sctrXiZeta), Uy(sctrXiZeta), Uz(sctrXiZeta)];    
                    xyz(:,i) = u + v;
                end
                plot_xyz_withMirroring
            end
        end
    end
    
    if plotEdge(3,1)
        % Plot surface where zeta = 0     
        zeta0Nodes = zeros(1,n*m);
        counter = 1;
        for j = 1:m
            for i = 1:n
                zeta0Nodes(counter) = n*(j-1) + i;
                counter = counter + 1;
            end
        end
        [XiEtaMesh, ~, ~, elRangeXi, elRangeEta] = generateIGA2DMesh(Xi, Eta, p, q, n, m);
        
     
        XYZ = zeros(length(Xi_values), length(Eta_values),3);
        W = zeros(length(Xi_values), length(Eta_values));
        for i = 1:length(Xi_values)
            for j = 1:length(Eta_values)
                xi = Xi_values(i);
                eta = Eta_values(j);
 
                [R_fun, ~, ~] = NURBS2DBasis(xi, eta, p, q, Xi, Eta, weights(zeta0Nodes));

                idx_xi = findElementInMesh(xi, elRangeXi);
                idx_eta = findElementInMesh(eta, elRangeEta);

                e = noElemsXi*(idx_eta-1) + idx_xi;

                sctrXiEta = zeta0Nodes(XiEtaMesh(e,:));          %  element scatter vector
                v = R_fun*controlPts(sctrXiEta,:);
                u = R_fun*[Ux(sctrXiEta), Uy(sctrXiEta), Uz(sctrXiEta)];
                
                switch scalarfield
                    case 1
                        W(i,j) = signFactor*u(1);
                    case 2
                        W(i,j) = signFactor*u(2);
                    case 3
                        W(i,j) = signFactor*u(3);
                    case 4
                        W(i,j) = norm(u);
                    case 5
                        x = v(1);
                        y = v(2);
                        u_x = signFactor*u(1);
                        u_y = signFactor*u(2);
                        theta = atan2(y,x);
                        W(i,j) = u_y*sin(theta) + u_x*cos(theta);
                    case 6
                        x = v(1);
                        y = v(2);
                        z = v(3);
                        radius = sqrt(x^2+y^2+z^2);
                        phi = atan2(y,x);
                        theta = acos(z/radius);

                        u_x = signFactor*u(1);
                        u_y = signFactor*u(2);
                        u_z = signFactor*u(3);

                        W(i,j) = sin(theta)*cos(phi)*u_x + sin(theta)*sin(phi)*u_y + cos(theta)*u_z;
                end
                XYZ(i,j,:) = u + v;
            end
        end
        if minW > min(min(W))
            minW = min(min(W));
        end
        if maxW < max(max(W))
            maxW = max(max(W));
        end
        plotXYZwithMirroring
        
        if plotElementEdges
            % Plot element mesh where zeta = 0
            xyz = zeros(3, length(Eta_values));
            for i = 1:length(unique_xiKnots)
                xi = unique_xiKnots(i);
                for j = 1:length(Eta_values)
                    eta = Eta_values(j);
                    [R_fun, ~, ~] = NURBS2DBasis(xi, eta, p, q, Xi, Eta, weights(zeta0Nodes));

                    idx_xi = findElementInMesh(xi, elRangeXi);
                    idx_eta = findElementInMesh(eta, elRangeEta);

                    e = noElemsXi*(idx_eta-1) + idx_xi;

                    sctrXiEta = zeta0Nodes(XiEtaMesh(e,:));          %  element scatter vector
                    v = R_fun*controlPts(sctrXiEta,:);
                    u = R_fun*[Ux(sctrXiEta), Uy(sctrXiEta), Uz(sctrXiEta)];
                    xyz(:,j) = u + v;
                end
                plot_xyz_withMirroring
            end

            xyz = zeros(3, length(Xi_values));
            for j = 1:length(unique_etaKnots)
                eta = unique_etaKnots(j);
                for i = 1:length(Xi_values)
                    xi = Xi_values(i);
                    [R_fun, ~, ~] = NURBS2DBasis(xi, eta, p, q, Xi, Eta, weights(zeta0Nodes));

                    idx_xi = findElementInMesh(xi, elRangeXi);
                    idx_eta = findElementInMesh(eta, elRangeEta);

                    e = noElemsXi*(idx_eta-1) + idx_xi;

                    sctrXiEta = zeta0Nodes(XiEtaMesh(e,:));          %  element scatter vector
                    v = R_fun*controlPts(sctrXiEta,:);
                    u = R_fun*[Ux(sctrXiEta), Uy(sctrXiEta), Uz(sctrXiEta)];
                    xyz(:,i) = u + v;
                end
                plot_xyz_withMirroring
            end
        end
    end
    if plotEdge(3,2)
        %% Plot surface where zeta = 1
        zeta1Nodes = zeros(1,n*m);
        counter = 1;
        for j = 1:m
            for i = 1:n
                zeta1Nodes(counter) = (m*n)*(l-1) + n*(j-1) + i;
                counter = counter + 1;
            end
        end

        [XiEtaMesh, ~, ~, elRangeXi, elRangeEta] = generateIGA2DMesh(Xi, Eta, p, q, n, m);
        
        XYZ = zeros(length(Xi_values), length(Eta_values),3);
        W = zeros(length(Xi_values), length(Eta_values));
        for i = 1:length(Xi_values)
            for j = 1:length(Eta_values)
                xi = Xi_values(i);
                eta = Eta_values(j);
                [R_fun, ~, ~] = NURBS2DBasis(xi, eta, p, q, Xi, Eta, weights(zeta1Nodes));

                idx_xi = findElementInMesh(xi, elRangeXi);
                idx_eta = findElementInMesh(eta, elRangeEta);

                e = noElemsXi*(idx_eta-1) + idx_xi;

                sctrXiEta = zeta1Nodes(XiEtaMesh(e,:));          %  element scatter vector
                v = R_fun*controlPts(sctrXiEta,:);
                u = R_fun*[Ux(sctrXiEta), Uy(sctrXiEta), Uz(sctrXiEta)];
                switch scalarfield
                    case 1
                        W(i,j) = signFactor*u(1);
                    case 2
                        W(i,j) = signFactor*u(2);
                    case 3
                        W(i,j) = signFactor*u(3);
                    case 4
                        W(i,j) = norm(u);
                    case 5
                        x = v(1);
                        y = v(2);
                        u_x = signFactor*u(1);
                        u_y = signFactor*u(2);
                        theta = atan2(y,x);
                        W(i,j) = u_y*sin(theta) + u_x*cos(theta);
                    case 6
                        x = v(1);
                        y = v(2);
                        z = v(3);
                        radius = sqrt(x^2+y^2+z^2);
                        phi = atan2(y,x);
                        theta = acos(z/radius);

                        u_x = signFactor*u(1);
                        u_y = signFactor*u(2);
                        u_z = signFactor*u(3);

                        W(i,j) = sin(theta)*cos(phi)*u_x + sin(theta)*sin(phi)*u_y + cos(theta)*u_z;
                end
                XYZ(i,j,:) = u + v;
                Y(i,j) = u(2) + v(2);
                Z(i,j) = u(3) + v(3);
            end
        end
        if minW > min(min(W))
            minW = min(min(W));
        end
        if maxW < max(max(W))
            maxW = max(max(W));
        end
        plotXYZwithMirroring

        if plotElementEdges
            % Plot element mesh where zeta = 1
            xyz = zeros(3, length(Eta_values));
            for i = 1:length(unique_xiKnots)
                xi = unique_xiKnots(i);
                for j = 1:length(Eta_values)
                    eta = Eta_values(j);
                    [R_fun, ~, ~] = NURBS2DBasis(xi, eta, p, q, Xi, Eta, weights(zeta1Nodes));

                    idx_xi = findElementInMesh(xi, elRangeXi);
                    idx_eta = findElementInMesh(eta, elRangeEta);

                    e = noElemsXi*(idx_eta-1) + idx_xi;

                    sctrXiEta = zeta1Nodes(XiEtaMesh(e,:));          %  element scatter vector
                    v = R_fun*controlPts(sctrXiEta,:);
                    u = R_fun*[Ux(sctrXiEta), Uy(sctrXiEta), Uz(sctrXiEta)];
                    xyz(:,j) = u + v;
                end
                plot_xyz_withMirroring
            end

            xyz = zeros(3, length(Xi_values));
            for j = 1:length(unique_etaKnots)
                eta = unique_etaKnots(j);
                for i = 1:length(Xi_values)
                    xi = Xi_values(i);
                    [R_fun, ~, ~] = NURBS2DBasis(xi, eta, p, q, Xi, Eta, weights(zeta1Nodes));

                    idx_xi = findElementInMesh(xi, elRangeXi);
                    idx_eta = findElementInMesh(eta, elRangeEta);

                    e = noElemsXi*(idx_eta-1) + idx_xi;

                    sctrXiEta = zeta1Nodes(XiEtaMesh(e,:));          %  element scatter vector
                    v = R_fun*controlPts(sctrXiEta,:);
                    u = R_fun*[Ux(sctrXiEta), Uy(sctrXiEta), Uz(sctrXiEta)];
                    xyz(:,i) = u + v;
                end
                plot_xyz_withMirroring
            end
        end
    end
    if abs(minW) > abs(maxW)
        maxAbsW = abs(minW);
    else
        maxAbsW = abs(maxW);
    end
    if minW > 0
        caxis([0 maxAbsW])
        customColorMap = customColorMap(((end+1)/2):end,:);
    elseif maxW < 0
        caxis([-maxAbsW 0])
        customColorMap = customColorMap(1:((end+1)/2),:);
    else
        caxis([-maxAbsW maxAbsW])
    end
    colormap(customColorMap)
    h = gcf;
    
elseif strcmp(nurbs.type, '2Dsurface')
    maxW = -inf;
    minW = inf;
    p = nurbs.degree(1);
    q = nurbs.degree(2);
    
    n = nurbs.number(1);
    m = nurbs.number(2);
    
    Xi = nurbs.knots{1};
    Eta = nurbs.knots{2};
    
    if iscell(results)
        Ux = results{1};
        Uy = results{2};
        U = [Ux, Uy];
    else
        U = results;
    end
    
    % To reconstruct any interpolatory points, any p repeated knot should be
    % included in the abscissa values to plot
    Xi_values = insertUniform3(nurbs.knots{1}, resolution(1), nurbs.degree(1));
    Eta_values = insertUniform3(nurbs.knots{2}, resolution(2), nurbs.degree(2));


    X = zeros(length(Xi_values), length(Eta_values));
    Y = zeros(length(Xi_values), length(Eta_values));
    W = zeros(length(Xi_values), length(Eta_values));
    
    for i = 1:length(Xi_values)
        disp([num2str(i) ' of ' num2str(length(Xi_values)) ' completed!'])
        for j = 1:length(Eta_values)
            xi = Xi_values(i);
            eta = Eta_values(j);
            [u, v] = numericalSolEval_2D_final(xi, eta, p, q, Xi, Eta, weights, controlPts, U, 0);
            W(i,j) = u;
%             XYZ(i,j,:) = u + v;
%             X(i,j) = u(2) + v(2);
%             Y(i,j) = u(3) + v(3);
            X(i,j) = v(1);
            Y(i,j) = v(2);
        end
    end
%     surf(X,Y,real(W), 'EdgeColor','none', 'LineStyle','none')
%     view(0,90)
%     hold on
    
%         plotXYZwithMirroring
    unique_xiKnots = unique(nurbs.knots{1});
    unique_etaKnots = unique(nurbs.knots{2});
    
    if plotElementEdges
        xy = zeros(3, length(Eta_values));
        for i = 1:length(unique_xiKnots)
            xi = unique_xiKnots(i);
            for j = 1:length(Eta_values)
                eta = Eta_values(j);
                [u, v] = numericalSolEval_2D_final(xi, eta, p, q, Xi, Eta, weights, controlPts, U, 0);
                xy(1:2,j) = v;
                xy(3,j) = u;
    %             xyz(:,j) = u + v;
            end
            plot3(xy(1,:),xy(2,:),real(xy(3,:)),'black')
        end

        xy = zeros(3, length(Xi_values));
        for j = 1:length(unique_etaKnots)
            eta = unique_etaKnots(j);
            for i = 1:length(Xi_values)
                xi = Xi_values(i);
                [u, v] = numericalSolEval_2D_final(xi, eta, p, q, Xi, Eta, weights, controlPts, U, 0);
                xy(1:2,i) = v;
                xy(3,i) = u;
    %             xyz(:,j) = u + v;
            end
            plot3(xy(1,:),xy(2,:),real(xy(3,:)),'black')
        end
    end
%     end

%     if abs(minW) > abs(maxW)
%         maxAbsW = abs(minW);
%     else
%         maxAbsW = abs(maxW);
%     end
%     if minW > 0
%         caxis([0 maxAbsW])
%         customColorMap = customColorMap(((end+1)/2):end,:);
%     elseif maxW < 0
%         caxis([-maxAbsW 0])
%         customColorMap = customColorMap(1:((end+1)/2),:);
%     else
%         caxis([-maxAbsW maxAbsW])
%     end
%     colormap(customColorMap)
    h = gcf;    
    
elseif strcmp(nurbs.type, '2Dcurve')
    Xi_values = insertUniform(nurbs.knots, resolution);

    C = zeros(length(Xi_values), 2);

    for j = 1:length(Xi_values)
        xi = Xi_values(j);
        C(j,:) = evaluateNURBS(nurbs, xi);
    end

    h = plot(C(:,1), C(:,2), 'color', mycolor);    
end






