function h = plotSpline(splineObj, resolution)

% Plots the surface of a NURBS volume
% colors = colormap('summer');
color = 1.5*[44 77 32]/255; % colors(end,:);
mycolor = [33 76 161]/255;
% color = [33 76 161]/255; %214CA1
lineColor = 'black';
if strcmp(splineObj.type, '3Dvolume')
    resU = resolution(1);
    resV = resolution(2);
    resW = resolution(3);

    unique_uKnots = unique(splineObj.knots{1});
    unique_vKnots = unique(splineObj.knots{2});
    unique_wKnots = unique(splineObj.knots{3});

    Xi = linspace(0,1,resU);
    Eta = linspace(0,1,resV);
    Zeta = linspace(0,1,resW);

    %% Plot surface where xi = 0
    X = zeros(resV, resW);
    Y = zeros(resV, resW);
    Z = zeros(resV, resW);
    for j = 1:resV
        for k = 1:resW
            xi = 0;
            eta = Eta(j);
            zeta = Zeta(k);
            v = evaluateSpline(splineObj, [xi eta zeta]);
            X(j,k) = v(1);
            Y(j,k) = v(2);
            Z(j,k) = v(3);
        end
    end
    surf(X,Y,Z, 'FaceColor', color,'EdgeColor','none','LineStyle','none')
    hold on

    % Plot element mesh where xi = 0
    x = zeros(resW,1);
    y = zeros(resW,1);
    z = zeros(resW,1);
    for j = 1:length(unique_vKnots)
        xi = 0;
        eta = unique_vKnots(j);
        for k = 1:resW
            zeta = Zeta(k);
            v = evaluateSpline(splineObj, [xi eta zeta]);
            x(k) = v(1);
            y(k) = v(2);
            z(k) = v(3);
        end
        if zeta == 0.25 || zeta == 0.5 || zeta == 0.75
            plot3(x,y,z,'color','red')
        else
            plot3(x,y,z,'color',lineColor)
        end  
    end

    x = zeros(resV,1);
    y = zeros(resV,1);
    z = zeros(resV,1);
    for k = 1:length(unique_wKnots)
        xi = 0;
        zeta = unique_wKnots(k);
        for j = 1:resV
            eta = Eta(j);
            v = evaluateSpline(splineObj, [xi eta zeta]);
            x(j) = v(1);
            y(j) = v(2);
            z(j) = v(3);
        end
        if zeta == 0.25 || zeta == 0.5 || zeta == 0.75
            plot3(x,y,z,'color','red')
        else
            plot3(x,y,z,'color',lineColor)
        end  
    end


    %% Plot surface where xi = 1
    X = zeros(resV, resW);
    Y = zeros(resV, resW);
    Z = zeros(resV, resW);
    for j = 1:resV
        for k = 1:resW
            xi = 1;
            eta = Eta(j);
            zeta = Zeta(k);
            v = evaluateSpline(splineObj, [xi eta zeta]);
            X(j,k) = v(1);
            Y(j,k) = v(2);
            Z(j,k) = v(3);
        end
    end
    surf(X,Y,Z, 'FaceColor', color,'EdgeColor','none','LineStyle','none')

    % Plot element mesh where xi = 1
    x = zeros(resW,1);
    y = zeros(resW,1);
    z = zeros(resW,1);
    for j = 1:length(unique_vKnots)
        xi = 1;
        eta = unique_vKnots(j);
        for k = 1:resW
            zeta = Zeta(k);
            v = evaluateSpline(splineObj, [xi eta zeta]);
            x(k) = v(1);
            y(k) = v(2);
            z(k) = v(3);
        end
        if zeta == 0.25 || zeta == 0.5 || zeta == 0.75
            plot3(x,y,z,'color','red')
        else
            plot3(x,y,z,'color',lineColor)
        end  
    end

    x = zeros(resV,1);
    y = zeros(resV,1);
    z = zeros(resV,1);
    for k = 1:length(unique_wKnots)
        xi = 1;
        zeta = unique_wKnots(k);
        for j = 1:resV
            eta = Eta(j);
            v = evaluateSpline(splineObj, [xi eta zeta]);
            x(j) = v(1);
            y(j) = v(2);
            z(j) = v(3);
        end
        if zeta == 0.25 || zeta == 0.5 || zeta == 0.75
            plot3(x,y,z,'color','red')
        else
            plot3(x,y,z,'color',lineColor)
        end  
    end

    %% Plot surface where eta = 0
    for i = 1:resU
        for k = 1:resW
            xi = Xi(i);
            eta = 0;
            zeta = Zeta(k);
            v = evaluateSpline(splineObj, [xi eta zeta]);
            X(i,k) = v(1);
            Y(i,k) = v(2);
            Z(i,k) = v(3);
        end
    end
    surf(X,Y,Z, 'FaceColor', color,'EdgeColor','none','LineStyle','none')

    % Plot element mesh where eta = 0
    x = zeros(resW,1);
    y = zeros(resW,1);
    z = zeros(resW,1);
    for i = 1:length(unique_uKnots)
        xi = unique_uKnots(i);
        eta = 0;
        for k = 1:resW
            zeta = Zeta(k);
            v = evaluateSpline(splineObj, [xi eta zeta]);
            x(k) = v(1);
            y(k) = v(2);
            z(k) = v(3);
        end
        if zeta == 0.25 || zeta == 0.5 || zeta == 0.75
            plot3(x,y,z,'color','red')
        else
            plot3(x,y,z,'color',lineColor)
        end  
    end

    x = zeros(resU,1);
    y = zeros(resU,1);
    z = zeros(resU,1);
    for k = 1:length(unique_wKnots)
        eta = 0;
        zeta = unique_wKnots(k);
        for i = 1:resU
            xi = Xi(i);
            v = evaluateSpline(splineObj, [xi eta zeta]);
            x(i) = v(1);
            y(i) = v(2);
            z(i) = v(3);
        end
        if zeta == 0.25 || zeta == 0.5 || zeta == 0.75
            plot3(x,y,z,'color','red')
        else
            plot3(x,y,z,'color',lineColor)
        end  
    end

    %% Plot surface where eta = 1
    X = zeros(resU, resW);
    Y = zeros(resU, resW);
    Z = zeros(resU, resW);
    for i = 1:resU
        for k = 1:resW
            xi = Xi(i);
            eta = 1;
            zeta = Zeta(k);
            v = evaluateSpline(splineObj, [xi eta zeta]);
            X(i,k) = v(1);
            Y(i,k) = v(2);
            Z(i,k) = v(3);
        end
    end
    surf(X,Y,Z, 'FaceColor', color,'EdgeColor','none','LineStyle','none')

    % Plot element mesh where eta = 1
    x = zeros(resW,1);
    y = zeros(resW,1);
    z = zeros(resW,1);
    for i = 1:length(unique_uKnots)
        xi = unique_uKnots(i);
        eta = 1;
        for k = 1:resW
            zeta = Zeta(k);
            v = evaluateSpline(splineObj, [xi eta zeta]);
            x(k) = v(1);
            y(k) = v(2);
            z(k) = v(3);
        end
        if zeta == 0.25 || zeta == 0.5 || zeta == 0.75
            plot3(x,y,z,'color','red')
        else
            plot3(x,y,z,'color',lineColor)
        end  
    end

    x = zeros(resU,1);
    y = zeros(resU,1);
    z = zeros(resU,1);
    for k = 1:length(unique_wKnots)
        eta = 1;
        zeta = unique_wKnots(k);
        for i = 1:resU
            xi = Xi(i);
            v = evaluateSpline(splineObj, [xi eta zeta]);
            x(i) = v(1);
            y(i) = v(2);
            z(i) = v(3);
        end
        if zeta == 0.25 || zeta == 0.5 || zeta == 0.75
            plot3(x,y,z,'color','red')
        else
            plot3(x,y,z,'color',lineColor)
        end  
    end

    %% Plot surface where zeta = 0
    X = zeros(resU, resV);
    Y = zeros(resU, resV);
    Z = zeros(resU, resV);
    for i = 1:resU
        for j = 1:resV
            xi = Xi(i);
            eta = Eta(j);
            zeta = 0;
            v = evaluateSpline(splineObj, [xi eta zeta]);
            X(i,j) = v(1);
            Y(i,j) = v(2);
            Z(i,j) = v(3);
        end
    end
    surf(X,Y,Z, 'FaceColor', color,'EdgeColor','none','LineStyle','none')

    % Plot element mesh where zeta = 0
    x = zeros(resV,1);
    y = zeros(resV,1);
    z = zeros(resV,1);
    for i = 1:length(unique_uKnots)
        xi = unique_uKnots(i);
        zeta = 0;
        for j = 1:resV
            eta = Eta(j);
            v = evaluateSpline(splineObj, [xi eta zeta]);
            x(j) = v(1);
            y(j) = v(2);
            z(j) = v(3);
        end
        if zeta == 0.25 || zeta == 0.5 || zeta == 0.75
            plot3(x,y,z,'color','red')
        else
            plot3(x,y,z,'color',lineColor)
        end  
    end

    x = zeros(resU,1);
    y = zeros(resU,1);
    z = zeros(resU,1);
    for j = 1:length(unique_vKnots)
        eta = unique_vKnots(j);
        zeta = 0;
        for i = 1:resU
            xi = Xi(i);
            v = evaluateSpline(splineObj, [xi eta zeta]);
            x(i) = v(1);
            y(i) = v(2);
            z(i) = v(3);
        end
        if zeta == 0.25 || zeta == 0.5 || zeta == 0.75
            plot3(x,y,z,'color','red')
        else
            plot3(x,y,z,'color',lineColor)
        end  
    end

    %% Plot surface where zeta = 1
    X = zeros(resU, resV);
    Y = zeros(resU, resV);
    Z = zeros(resU, resV);
    for i = 1:resU
        for j = 1:resV
            xi = Xi(i);
            eta = Eta(j);
            zeta = 1;
            v = evaluateSpline(splineObj, [xi eta zeta]);
            X(i,j) = v(1);
            Y(i,j) = v(2);
            Z(i,j) = v(3);
        end
    end
    surf(X,Y,Z, 'FaceColor', color,'EdgeColor','none','LineStyle','none');

    %%%%%%%%%%%%% HER?
    % Plot element mesh where zeta = 1
    x = zeros(resV,1);
    y = zeros(resV,1);
    z = zeros(resV,1);
    for i = 1:length(unique_uKnots)
        xi = unique_uKnots(i);
        zeta = 1;
        for j = 1:resV
            eta = Eta(j);
            v = evaluateSpline(splineObj, [xi eta zeta]);
            x(j) = v(1);
            y(j) = v(2);
            z(j) = v(3);
        end
        if zeta == 0.25 || zeta == 0.5 || zeta == 0.75
            plot3(x,y,z,'color','red')
        else
            plot3(x,y,z,'color',lineColor)
        end  
    end

    x = zeros(resU,1);
    y = zeros(resU,1);
    z = zeros(resU,1);
    for j = 1:length(unique_vKnots)
        eta = unique_vKnots(j);
        zeta = 1;
        for i = 1:resU
            xi = Xi(i);
            v = evaluateSpline(splineObj, [xi eta zeta]);
            x(i) = v(1);
            y(i) = v(2);
            z(i) = v(3);
        end
        if zeta == 0.25 || zeta == 0.5 || zeta == 0.75
            plot3(x,y,z,'color','red')
        else
            plot3(x,y,z,'color',lineColor)
        end  
    end

    daspect([1 1 1]);
    view(3)
    camlight;
    lighting phong;
    set(gca, 'Color', 'none');
    camlight('left')
    axis equal
    h = gcf;
elseif strcmp(splineObj.type, '2Dcurve')
    Xi = splineObj.knots;
    xi_array = linspace(Xi(1),Xi(end),resolution);

    C = zeros(length(xi_array), 2);

    for j = 1:length(xi_array)
        xi = xi_array(j);
        C(j,:) = evaluateSpline(splineObj, xi);
    end

    h = plot(C(:,1), C(:,2), 'color', mycolor);
    
end





