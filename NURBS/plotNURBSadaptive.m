function [minC, maxC] = plotNURBSadaptive(nurbs, resolution, plotElementEdges, U, varCol, Eta_e, Xi_e, minC, maxC)
if nargin < 8
    minC = inf;
    maxC = -inf;
end
p_xi = nurbs.degree(1);
p_eta = nurbs.degree(2);
Xi = nurbs.knots{1};
Eta = nurbs.knots{2};
weights = varCol.weights;
controlPts = varCol.controlPts;
colorFun = varCol.colorFun;
plotRealPart = false;
if plotRealPart
    maxValue = 2;
else
    maxValue = varCol.k/(4*pi)+10*eps;
end
x = varCol.x;
if nargin < 6
    hold on
    axis equal
    axis off
    set(gca, 'Color', 'none');
    view(-90,-30)
    colorbar
    drawnow
    hold on
    uniqueXi = unique(Xi);
    uniqueEta = unique(Eta);
    
    n_div_xi = length(uniqueXi)-1;
    n_div_eta = length(uniqueEta)-1;
    for i_eta = 1:n_div_eta
        Eta_e_sub = uniqueEta(i_eta:i_eta+1);
        for i_xi = 1:n_div_xi
            Xi_e_sub = uniqueXi(i_xi:i_xi+1);
            [minC, maxC] = plotNURBSadaptive(nurbs, resolution, plotElementEdges, U, varCol, Eta_e_sub, Xi_e_sub, minC, maxC);
        end
    end
    Xi_values = insertUniform3(nurbs.knots{1}, 10*resolution(1), nurbs.degree(1));
    Eta_values = insertUniform3(nurbs.knots{2}, 10*resolution(2), nurbs.degree(2));
    if plotElementEdges        
        X = zeros(length(Xi_values),1);
        Y = zeros(length(Xi_values),1);
        Z = zeros(length(Xi_values),1);
        X2 = zeros(length(Xi_values),1);
        Y2 = zeros(length(Xi_values),1);
        Z2 = zeros(length(Xi_values),1);

        for j = 1:length(uniqueEta)
            eta = uniqueEta(j);
            for i = 1:length(Xi_values)
                xi = Xi_values(i);
                v = evaluateNURBS(nurbs, [xi eta]);
                X(i) = v(1);
                Y(i) = v(2);
                Z(i) = v(3);
                
                
                xmy = x - v.';
                nx = x'/norm(x');
                ny = v/norm(v);
                r = norm(xmy);
                if plotRealPart
                    f_eval = real(colorFun(xmy,r,nx,ny));
                else
                    f_eval = imag(colorFun(xmy,r,nx,ny));
                end
                    
                X(j) = v(1);
                Y(j) = v(2);
                Z(j) = v(3);
                if f_eval > maxValue
                    f_eval = NaN;
                end
                X2(i) = v(1) + v(1)/norm(v)*f_eval/maxValue*norm(v)/3;
                Y2(i) = v(2) + v(2)/norm(v)*f_eval/maxValue*norm(v)/3;
                Z2(i) = v(3) + v(3)/norm(v)*f_eval/maxValue*norm(v)/3;
            end
%             plot3(X,Y,Z,'color','black')
            plot3(X2,Y2,Z2,'color','black')
        end

        X = zeros(length(Eta_values),1);
        Y = zeros(length(Eta_values),1);
        Z = zeros(length(Eta_values),1);
        X2 = zeros(length(Eta_values),1);
        Y2 = zeros(length(Eta_values),1);
        Z2 = zeros(length(Eta_values),1);
        for i = 1:length(uniqueXi)
            xi = uniqueXi(i);
            for j = 1:length(Eta_values)
                eta = Eta_values(j);
                v = evaluateNURBS(nurbs, [xi eta]);
                X(j) = v(1);
                Y(j) = v(2);
                Z(j) = v(3);
                
                
                xmy = x - v.';
                nx = x'/norm(x');
                ny = v/norm(v);
                r = norm(xmy);
                if plotRealPart
                    f_eval = real(colorFun(xmy,r,nx,ny));
                else
                    f_eval = imag(colorFun(xmy,r,nx,ny));
                end
                X(j) = v(1);
                Y(j) = v(2);
                Z(j) = v(3);
                if f_eval > maxValue
                    f_eval = NaN;
                end
                X2(j) = v(1) + v(1)/norm(v)*f_eval/maxValue*norm(v)/3;
                Y2(j) = v(2) + v(2)/norm(v)*f_eval/maxValue*norm(v)/3;
                Z2(j) = v(3) + v(3)/norm(v)*f_eval/maxValue*norm(v)/3;
            end
%             plot3(X,Y,Z,'color','black')
            plot3(X2,Y2,Z2,'color','black')
        end
    end
    hold off
else
    
    h_min = varCol.h_min;
    x_1 = evaluateNURBS(nurbs, [Xi_e(1)+eps, Eta_e(1)+eps]).';
    x_2 = evaluateNURBS(nurbs, [Xi_e(2)-eps, Eta_e(1)+eps]).';
    x_3 = evaluateNURBS(nurbs, [Xi_e(2)-eps, Eta_e(2)-eps]).';
    x_4 = evaluateNURBS(nurbs, [Xi_e(1)+eps, Eta_e(2)-eps]).';
    % x_5 = evaluateNURBS(nurbs, [mean(Xi_e),  mean(Eta_e)]).';

    l = min([norm(x-x_1), norm(x-x_2), norm(x-x_3), norm(x-x_4)]);
    h_1 = norm(x_1-x_3);
    h_2 = norm(x_2-x_4);
    h = max(h_1,h_2);

    coeff2 = 2;
    coeff3 = 0.05;


    xi_h = max(norm(x_1-x_2),norm(x_3-x_4));
    eta_h = max(norm(x_2-x_3),norm(x_1-x_4));
    
%     keyboard
    if h/l > coeff2 && h/h_min > coeff3
    %     if h/h_min > coeff3 %|| ((xi_h/eta_h > 2 || eta_h/xi_h > 2) && h/h_min > 0.5*coeff3)
        if xi_h/eta_h > 2
            n_div_xi = round(xi_h/eta_h);
            n_div_eta = 1;
        elseif eta_h/xi_h > 2
            n_div_xi = 1;
            n_div_eta = round(eta_h/xi_h);
        else
            n_div_xi = 2;
            n_div_eta = 2;
        end

        Xi_e_arr  = linspace(Xi_e(1),Xi_e(2),n_div_xi+1);
        Eta_e_arr = linspace(Eta_e(1),Eta_e(2),n_div_eta+1);

        for i_eta = 1:n_div_eta
            Eta_e_sub = Eta_e_arr(i_eta:i_eta+1);
            for i_xi = 1:n_div_xi
                Xi_e_sub = Xi_e_arr(i_xi:i_xi+1);
                [minC, maxC] = plotNURBSadaptive(nurbs, resolution, plotElementEdges, U, varCol, Eta_e_sub, Xi_e_sub, minC, maxC);
            end
        end
    else
        Xi_values = linspace(Xi_e(1), Xi_e(2), resolution(1));
        Eta_values = linspace(Eta_e(1), Eta_e(2), resolution(2));
        X = zeros(length(Xi_values), length(Eta_values));
        Y = zeros(length(Xi_values), length(Eta_values));
        Z = zeros(length(Xi_values), length(Eta_values));
        X2 = zeros(length(Xi_values), length(Eta_values));
        Y2 = zeros(length(Xi_values), length(Eta_values));
        Z2 = zeros(length(Xi_values), length(Eta_values));
        C = zeros(length(Xi_values), length(Eta_values));
        parfor i = 1:length(Xi_values)
            X_temp = zeros(1,length(Eta_values));
            Y_temp = zeros(1,length(Eta_values));
            Z_temp = zeros(1,length(Eta_values));
            X2_temp = zeros(1,length(Eta_values));
            Y2_temp = zeros(1,length(Eta_values));
            Z2_temp = zeros(1,length(Eta_values));
            C_temp = zeros(1,length(Eta_values));
            xi = Xi_values(i);
            for j = 1:length(Eta_values)
                eta = Eta_values(j);
                v = evaluateNURBS(nurbs, [xi eta]);
                xmy = x - v.';
                nx = x'/norm(x');
                ny = v/norm(v);
                r = norm(xmy);
                if plotRealPart
                    f_eval = real(colorFun(xmy,r,nx,ny));
                else
                    f_eval = imag(colorFun(xmy,r,nx,ny));
                end
                if f_eval > maxValue
                    f_eval = NaN;
                end
                C_temp(j) = f_eval;
                X_temp(j) = v(1);
                Y_temp(j) = v(2);
                Z_temp(j) = v(3);
                X2_temp(j) = v(1) + v(1)/norm(v)*f_eval/maxValue*norm(v)/3;
                Y2_temp(j) = v(2) + v(2)/norm(v)*f_eval/maxValue*norm(v)/3;
                Z2_temp(j) = v(3) + v(3)/norm(v)*f_eval/maxValue*norm(v)/3;
            end
            X(i,:) = X_temp;
            Y(i,:) = Y_temp;
            Z(i,:) = Z_temp;
            X2(i,:) = X2_temp;
            Y2(i,:) = Y2_temp;
            Z2(i,:) = Z2_temp;
            C(i,:) = C_temp;
        end
        if plotElementEdges
            surf(X,Y,Z, 'FaceColor', 'white', 'EdgeColor','none','LineStyle','none','FaceAlpha',0.2)
            surf(X2,Y2,Z2, C, 'EdgeColor','none','LineStyle','none')
        else
            surf(X,Y,Z, C)
            surf(X2,Y2,Z2, C)
        end
        if minC > min(min(C))
            minC = min(min(C));
        end
        if maxC < max(max(C))
            maxC = max(max(C));
        end
        caxis([minC,maxC])
    end
end