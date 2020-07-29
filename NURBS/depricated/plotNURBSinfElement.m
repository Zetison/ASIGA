function h = plotNURBSinfElement(nurbs, resolution, plotElementEdges, color, alphaValue,c_xy,c_z, x_0)
error('Depricated')
if nargin < 3
    plotElementEdges = true;
end
if nargin < 4
    % Plots the surface of a NURBS volume
    % colors = colormap('summer');
    color = getColor(1); % colors(end,:);
    mycolor = [33 76 161]/255;
end
if nargin < 5
    alphaValue = 1;
end
% color = [33 76 161]/255; %214CA1
lineColor = 'black';
% A_2 maps the x-axis to the z-axis, the y-axis to the x-axis, and
% the z-axis to the y-axis
A_2 = [0 1 0;
      0 0 1;
      1 0 0];
% A_2 = eye(3);

unique_xiKnots = unique(nurbs.knots{1});
unique_etaKnots = unique(nurbs.knots{2});

% To reconstruct any interpolatory points, any p repeated knot should be
% included in the abscissa values to plot
Xi_values = insertUniform3(nurbs.knots{1}, resolution(1), nurbs.degree(1));
Eta_values = insertUniform3(nurbs.knots{2}, resolution(2), nurbs.degree(2));
%     %% Plot surface where zeta = 1
X = zeros(length(Xi_values), length(Eta_values));
Y = zeros(length(Xi_values), length(Eta_values));
Z = zeros(length(Xi_values), length(Eta_values));
for i = 1:length(Xi_values)
    for j = 1:length(Eta_values)
        xi = Xi_values(i);
        eta = Eta_values(j);
        v = evaluateNURBS(nurbs, [xi eta]);
        X(i,j) = v(1);
        Y(i,j) = v(2);
        Z(i,j) = v(3);
    end
end
if plotElementEdges
    if alphaValue < 1
        surf(X,Y,Z, 'FaceColor', color,'EdgeColor','none','LineStyle','none','FaceAlpha',alphaValue, 'FaceLighting', 'none')
    else
        surf(X,Y,Z, 'FaceColor', color,'EdgeColor','none','LineStyle','none','FaceAlpha',alphaValue)
    end
else
    surf(X,Y,Z, 'FaceColor', color)
end
hold on

if plotElementEdges        
    x = zeros(length(Xi_values),1);
    y = zeros(length(Xi_values),1);
    z = zeros(length(Xi_values),1);

    for j = 1:length(unique_etaKnots)
        eta = unique_etaKnots(j);
        for i = 1:length(Xi_values)
            xi = Xi_values(i);
            v = evaluateNURBS(nurbs, [xi eta]);
            x(i) = v(1);
            y(i) = v(2);
            z(i) = v(3);
        end
        plot3(x,y,z,'color',lineColor)
    end

    x = zeros(length(Eta_values),1);
    y = zeros(length(Eta_values),1);
    z = zeros(length(Eta_values),1);
    for i = 1:length(unique_xiKnots)
        xi = unique_xiKnots(i);
        for j = 1:length(Eta_values)
            eta = Eta_values(j);
            v = evaluateNURBS(nurbs, [xi eta]);
            x(j) = v(1);
            y(j) = v(2);
            z(j) = v(3);
        end
        plot3(x,y,z,'color',lineColor)
    end
end

Upsilon = sqrt(c_z^2-c_xy^2);
[R_a, ~, ~] = evaluateProlateCoords(0,0,c_z,Upsilon);
fact = 1/4;
r_arr = [linspace(R_a, (1+2*fact)*R_a, 100);
         linspace((1+2*fact)*R_a, (1+4*fact)*R_a, 100)];
xi_start = unique_xiKnots(round(length(unique_xiKnots)*6/16));
xi_end = unique_xiKnots(round(length(unique_xiKnots)*6/16)+1);
eta_start = unique_etaKnots(round(length(unique_etaKnots)*14/16)-3);
eta_end = unique_etaKnots(round(length(unique_etaKnots)*14/16)-2);

xi_arr = [xi_start xi_end];
eta_arr = [eta_start eta_end];

theta_arr = zeros(4,1);
phi_arr = zeros(4,1);
counter = 1;

for i = 1:2
    xi = xi_arr(i);
    for j = 1:2
        eta = eta_arr(j);
        X = evaluateNURBS(nurbs, [xi eta]);
        
        Xt = A_2*(X-x_0);
        [~, theta, phi] = evaluateProlateCoords(Xt(1),Xt(2),Xt(3),Upsilon);
        theta_arr(counter) = theta;
        phi_arr(counter) = phi;
        counter = counter + 1;
    end
end
        
for i = 1:4
    theta = theta_arr(i);
    phi = phi_arr(i);
    for j = 1:2
        xt = sqrt(r_arr(j,:).^2-Upsilon^2)*sin(theta)*cos(phi);
        yt = sqrt(r_arr(j,:).^2-Upsilon^2)*sin(theta)*sin(phi);
        zt = r_arr(j,:)*cos(theta);
        Xt = [xt; yt; zt];
        X = A_2\Xt + x_0;
        x = X(1,:);
        y = X(2,:);
        z = X(3,:);
        if j == 1
            plot3(x,y,z,'color',lineColor)
        else
            plot3(x,y,z,'--', 'color',lineColor)
        end
    end
end

Xi_values = linspace(xi_start, xi_end,resolution(1));
Eta_values = linspace(eta_start, eta_end,resolution(2));

X = zeros(length(Xi_values), length(Eta_values));
Y = zeros(length(Xi_values), length(Eta_values));
Z = zeros(length(Xi_values), length(Eta_values));

alphaValue = 0.8;

for n = [0,1,2,4]
    for i = 1:length(Xi_values)
        xi = Xi_values(i);
        for j = 1:length(Eta_values)
            eta = Eta_values(j);
            xx = evaluateNURBS(nurbs, [xi eta]);
            Xt = A_2*(xx-x_0);
            [~, theta, phi] = evaluateProlateCoords(Xt(1),Xt(2),Xt(3),Upsilon);
            r = (1+n*fact)*R_a+1e-3;
            xt = sqrt(r^2-Upsilon^2)*sin(theta)*cos(phi);
            yt = sqrt(r^2-Upsilon^2)*sin(theta)*sin(phi);
            zt = r*cos(theta);
            Xt = [xt; yt; zt];
            xx = A_2\Xt + x_0;
            x = xx(1);
            y = xx(2);
            z = xx(3);
            
            X(i,j) = x;
            Y(i,j) = y;
            Z(i,j) = z;
        end
    end
%     if n > 1
%         alphaValue = 0.6;
%     end
    if n == 4
        surf(X,Y,Z, 'FaceColor', [0 0 0.8],'EdgeColor','none','LineStyle','none','FaceAlpha',alphaValue)
    else
        surf(X,Y,Z, 'FaceColor', [0.8 0 0],'EdgeColor','none','LineStyle','none','FaceAlpha',alphaValue)
    end
%     keyboard
%     if n == 1
%     else
%         surf(Z,X,Y, 'FaceColor', [0.8 0 0],'EdgeColor',[84 84 84]/256,'FaceAlpha',alphaValue)
%     end      
    X = zeros(length(Xi_values),1);
    Y = zeros(length(Xi_values),1);
    Z = zeros(length(Xi_values),1);

    for j = 1:length(eta_arr)
        eta = eta_arr(j);
        for i = 1:length(Xi_values)
            xi = Xi_values(i);
            xx = evaluateNURBS(nurbs, [xi eta]);
            Xt = A_2*(xx-x_0);
            [~, theta, phi] = evaluateProlateCoords(Xt(1),Xt(2),Xt(3),Upsilon);
            r = (1+n*fact)*R_a+1e-3;
            xt = sqrt(r^2-Upsilon^2)*sin(theta)*cos(phi);
            yt = sqrt(r^2-Upsilon^2)*sin(theta)*sin(phi);
            zt = r*cos(theta);
            Xt = [xt; yt; zt];
            xx = A_2\Xt + x_0;
            x = xx(1);
            y = xx(2);
            z = xx(3);
            X(i) = x;
            Y(i) = y;
            Z(i) = z;
        end
        plot3(X,Y,Z,'color',lineColor)
    end

    X = zeros(length(Eta_values),1);
    Y = zeros(length(Eta_values),1);
    Z = zeros(length(Eta_values),1);
    for i = 1:length(xi_arr)
        xi = xi_arr(i);
        for j = 1:length(Eta_values)
            eta = Eta_values(j);
            xx = evaluateNURBS(nurbs, [xi eta]);
            Xt = A_2*(xx-x_0);
            [~, theta, phi] = evaluateProlateCoords(Xt(1),Xt(2),Xt(3),Upsilon);
            r = (1+n*fact)*R_a+1e-3;
            xt = sqrt(r^2-Upsilon^2)*sin(theta)*cos(phi);
            yt = sqrt(r^2-Upsilon^2)*sin(theta)*sin(phi);
            zt = r*cos(theta);
            Xt = [xt; yt; zt];
            xx = A_2\Xt + x_0;
            x = xx(1);
            y = xx(2);
            z = xx(3);
            X(j) = x;
            Y(j) = y;
            Z(j) = z;
        end
        plot3(X,Y,Z,'color',lineColor)
    end
        
end

h = gcf;

hold off






