function h = plotExample_InfiniteElements(nurbs, resolution, plotElementEdges)


color = [66 116 48]/255; % colors(end,:);
myColor = [33 76 161]/255;
water = [187 217 238]/255;
myRed = [0.5 0 0];
myGreen = [0 0.5 0];

lineColor = 'black';

unique_xiKnots = unique(nurbs.knots{1});
unique_etaKnots = unique(nurbs.knots{2});

% To reconstruct any interpolatory points, any p repeated knot should be
% included in the abscissa values to plot
Xi_values = insertUniform3(nurbs.knots{1}, resolution(1), nurbs.degree(1));
Eta_values = insertUniform3(nurbs.knots{2}(find(nurbs.knots{2} <= 0.5)), resolution(2), nurbs.degree(2));
Eta_values2 = insertUniform3(nurbs.knots{2}(find(nurbs.knots{2} >= 0.5)), resolution(2), nurbs.degree(2));

X = zeros(2*length(Xi_values)+2*length(Eta_values)-3,1);
Y = zeros(2*length(Xi_values)+2*length(Eta_values)-3,1);
X2 = zeros(2*length(Xi_values)+2*length(Eta_values2)-3,1);
Y2 = zeros(2*length(Xi_values)+2*length(Eta_values2)-3,1);
    
scaling = 1.05;
for i = 1:length(unique_xiKnots)
    xi = unique_xiKnots(i);
    v = evaluateNURBS(nurbs, [xi 1]);
    x = v(1);
    y = v(2);
    f = sqrt(nurbs.coeffs(1,1,end)^2-nurbs.coeffs(2,3,end)^2);
    R_a = evaluateEllipticCoords(x,y,f);
    theta = atan2(R_a*y,x*sqrt(R_a^2-f^2));
    r = linspace(R_a,scaling*R_a,resolution(2));
    plot(r*cos(theta), sqrt(r.^2-f^2)*sin(theta), 'color', lineColor)
    hold on
    r = linspace(scaling*R_a,1.1*scaling*R_a,resolution(2));
    plot(r*cos(theta), sqrt(r.^2-f^2)*sin(theta), '--', 'color', lineColor)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

% along eta = 0.5
for i = 2:length(Xi_values)
    xi = Xi_values(i);
    eta = 0.5;
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
counter = 1;
% along eta = 0.5
for i = 1:length(Xi_values)
    xi = Xi_values(i);
    eta = 0.5;
    v = evaluateNURBS(nurbs, [xi eta]);
    X2(counter) = v(1);
    Y2(counter) = v(2);
    counter = counter + 1;
end
% along xi = 1
for i = 2:length(Eta_values2)
    xi = 1;
    eta = Eta_values2(i);
    v = evaluateNURBS(nurbs, [xi eta]);
    X2(counter) = v(1);
    Y2(counter) = v(2);
    counter = counter + 1;
end
% reverse order of arrays
Xi_values = Xi_values(end:-1:1);
Eta_values2 = Eta_values2(end:-1:1);

% along eta = 1
for i = 2:length(Xi_values)
    xi = Xi_values(i);
    eta = 1;
    v = evaluateNURBS(nurbs, [xi eta]);
    X2(counter) = v(1);
    Y2(counter) = v(2);
    counter = counter + 1;
end
% along xi = 0
for i = 2:length(Eta_values2)
    xi = 0;
    eta = Eta_values2(i);
    v = evaluateNURBS(nurbs, [xi eta]);
    X2(counter) = v(1);
    Y2(counter) = v(2);
    counter = counter + 1;
end    

% reverse back order of arrays
Xi_values = Xi_values(end:-1:1);
Eta_values2 = Eta_values2(end:-1:1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if plotElementEdges
    fill(X,Y, color,'EdgeColor','none','LineStyle','none')
    hold on
    fill(X2,Y2, water,'EdgeColor','none','LineStyle','none')
else
    fill(X,Y, color)
end
hold on
Eta_values = insertUniform3(nurbs.knots{2}(find(nurbs.knots{2} >= 0.5)), resolution(2), nurbs.degree(2));
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
        if eta == 0.5
            plot(x,y,'color',myRed)
        elseif eta == 1
            plot(x,y,'color',myColor)
        else
            plot(x,y,'color',lineColor)
        end
    end
%     
%     % Plot element mesh for constant xi
%     x = zeros(length(Xi_values),1);
%     y = zeros(length(Xi_values),1);
% 
%     eta = 1;
%     for i = 1:length(Xi_values)
%         xi = Xi_values(i);
%         v = 0.9*scaling*evaluateNURBS(nurbs, [xi eta]);
%         x(i) = v(1);
%         y(i) = v(2);
%     end
%     plot(x,y, ':','color',lineColor)
end

