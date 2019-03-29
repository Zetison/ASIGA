function [h,maxC] = plotLagrange(nurbs, resolution, plotElementEdges, colorNew, alphaValue, U, varCol)
if nargin < 3
    plotElementEdges = true;
end
% Plots the surface of a NURBS volume
% colors = colormap('summer');
color = 1.5*[44 77 32]/255; % colors(end,:);
% color = [33 76 161]/255;
if nargin >= 4 && ~any(isnan(colorNew))
    color = colorNew;
end
if nargin > 5 && isfield(varCol,'colorFun')
    plotSolution = 1;
    colorFun = varCol.colorFun;
%     alphaValue = 1;
else
    plotSolution = 0;
    colorFun = @(v) NaN;
end
if nargin < 5
    alphaValue = 1;
end
if nargin < 7
    varCol.elementBasedSamples = false;
    if ~isfield(varCol,'plotAt')
        varCol.plotAt = [1 1;
                         1 1;
                         1 1];
    end
end
% keyboard
% plotSolution = 0;
% color = [33 76 161]/255; %214CA1
lineColor = 'black';
if strcmp(nurbs.type, '3Dvolume')
    plotAt = varCol.plotAt;
    
    % To reconstruct any interpolatory points, any p repeated knot should be
    % included in the abscissa values to plot
    xi = linspace(-1,1, resolution(1));
    eta = linspace(-1,1, resolution(2));
    zeta = linspace(-1,1, resolution(3));

    if plotAt(3,1)
        Nxi = nurbs.number(1);
        Neta = nurbs.number(2);
        Nzeta = nurbs.number(3);
        c = nurbs.coeffs;
        Bxi = zeros(Nxi,resolution(1));
        Beta = zeros(Neta,resolution(2));
        GLL = nurbs.GLL;
        for i = 1:Nxi
            Bxi(i,:) = lagrangePolynomials(xi,i,Nxi,GLL{1});
        end
        for j = 1:Neta
            Beta(j,:) = lagrangePolynomials(eta,j,Neta,GLL{2});
        end
        nxi = length(xi);
        neta = length(eta);

        X = (Bxi.'*reshape(c(1,:,:,1),Nxi,Neta))*Beta;
        Y = (Bxi.'*reshape(c(2,:,:,1),Nxi,Neta))*Beta;
        Z = (Bxi.'*reshape(c(3,:,:,1),Nxi,Neta))*Beta;
        
        %% Plot surface where zeta = 0
        
        colorXi = 'red';
        colorEta = 'red';
        colorZeta = 'red';
        markerColor = 'black';
        markerEdgeColor = markerColor;
        
%         plot3(reshape(x(1,:,:,1),(Nxi+1)*(Neta+1),1),reshape(x(2,:,:,1),(Nxi+1)*(Neta+1),1),reshape(x(3,:,:,1),(Nxi+1)*(Neta+1),1),'o','MarkerFaceColor',markerColor, 'MarkerEdgeColor', markerEdgeColor)
%         surf(X,Y,Z, C, 'EdgeColor','none','LineStyle','none')
        hold on
%         surf(X,Y,Z, 'EdgeColor','none','LineStyle','none')
%         axis equal
        if plotSolution
            if ~any(isnan(U))
                u_h = (Bxi.'*U(:,:,1))*Beta;
%                 C = reshape(log10(abs(norm2([reshape(X,nxi*neta,1),reshape(Y,nxi*neta,1),reshape(Z,nxi*neta,1)])-1)),nxi,neta);
                if plotElementEdges
                    surf(X,Y,Z, real(u_h), 'EdgeColor','none','LineStyle','none')
                else
                    surf(X,Y,Z, real(u_h))
                end
            else
                C = reshape(log10(abs(norm2([reshape(X,nxi*neta,1),reshape(Y,nxi*neta,1),reshape(Z,nxi*neta,1)])-1)),nxi,neta);
                if plotElementEdges
                    surf(X,Y,Z, C, 'EdgeColor','none','LineStyle','none')
                else
                    surf(X,Y,Z, C)
                end
                colorbar
                caxis([-17,0])
                colorbar off
                [Xs,Ys,Zs] = sphere(100);
                for j = 1:Neta
                    for i = 1:Nxi
                        surf(1e-2*Xs+c(1,i,j,1),1e-2*Ys+c(2,i,j,1),1e-2*Zs+c(3,i,j,1), 'FaceColor', 'black','EdgeColor','none','LineStyle','none')
                    end
                end
            end
        else
            if true
                if alphaValue < 1
                    surf(X,Y,Z, 'FaceColor', color,'EdgeColor','none','LineStyle','none','FaceAlpha',1, 'FaceLighting', 'none')
                else
                    surf(X,Y,Z, 'FaceColor', color,'EdgeColor','none','LineStyle','none','FaceAlpha',1)
                end
            else
                surf(X,Y,Z, 'FaceColor', color)
            end
        end
        if plotElementEdges
            plot3(X(:,1),Y(:,1),Z(:,1),'color',lineColor)
            plot3(X(:,end),Y(:,end),Z(:,end),'color',lineColor)

            plot3(X(1,:),Y(1,:),Z(1,:),'color',lineColor)
            plot3(X(end,:),Y(end,:),Z(end,:),'color',lineColor)
        end
        
    end
    if plotAt(3,2)
        Nxi = nurbs.number(1);
        Neta = nurbs.number(2);
        Nzeta = nurbs.number(3);
        c = nurbs.coeffs;
        Bxi = zeros(Nxi,resolution(1));
        Beta = zeros(Neta,resolution(2));
        GLL = nurbs.GLL;
        for i = 1:Nxi
            Bxi(i,:) = lagrangePolynomials(xi,i,Nxi,GLL{1});
        end
        for j = 1:Neta
            Beta(j,:) = lagrangePolynomials(eta,j,Neta,GLL{2});
        end
        nxi = length(xi);
        neta = length(eta);

        X = (Bxi.'*reshape(c(1,:,:,end),Nxi,Neta))*Beta;
        Y = (Bxi.'*reshape(c(2,:,:,end),Nxi,Neta))*Beta;
        Z = (Bxi.'*reshape(c(3,:,:,end),Nxi,Neta))*Beta;
        
        %% Plot surface where zeta = 0
        
        colorXi = 'red';
        colorEta = 'red';
        colorZeta = 'red';
        markerColor = 'black';
        markerEdgeColor = markerColor;
        
%         plot3(reshape(x(1,:,:,1),(Nxi+1)*(Neta+1),1),reshape(x(2,:,:,1),(Nxi+1)*(Neta+1),1),reshape(x(3,:,:,1),(Nxi+1)*(Neta+1),1),'o','MarkerFaceColor',markerColor, 'MarkerEdgeColor', markerEdgeColor)
%         surf(X,Y,Z, C, 'EdgeColor','none','LineStyle','none')
        hold on
%         surf(X,Y,Z, 'EdgeColor','none','LineStyle','none')
%         axis equal
        if plotSolution
            if ~any(isnan(U))
                u_h = (Bxi.'*U(:,:,1))*Beta;
%                 C = reshape(log10(abs(norm2([reshape(X,nxi*neta,1),reshape(Y,nxi*neta,1),reshape(Z,nxi*neta,1)])-1)),nxi,neta);
                if plotElementEdges
                    surf(X,Y,Z, real(u_h), 'EdgeColor','none','LineStyle','none')
                else
                    surf(X,Y,Z, real(u_h))
                end
            else
                C = reshape(log10(abs(norm2([reshape(X,nxi*neta,1),reshape(Y,nxi*neta,1),reshape(Z,nxi*neta,1)])-1)),nxi,neta);
                if plotElementEdges
                    surf(X,Y,Z, C, 'EdgeColor','none','LineStyle','none')
                else
                    surf(X,Y,Z, C)
                end
                colorbar
                caxis([-17,0])
                colorbar off
                [Xs,Ys,Zs] = sphere(100);
                for j = 1:Neta
                    for i = 1:Nxi
                        surf(1e-2*Xs+c(1,i,j,1),1e-2*Ys+c(2,i,j,1),1e-2*Zs+c(3,i,j,1), 'FaceColor', 'black','EdgeColor','none','LineStyle','none')
                    end
                end
            end
        else
            if true
                if alphaValue < 1
                    surf(X,Y,Z, 'FaceColor', color,'EdgeColor','none','LineStyle','none','FaceAlpha',alphaValue, 'FaceLighting', 'none')
                else
                    surf(X,Y,Z, 'FaceColor', color,'EdgeColor','none','LineStyle','none','FaceAlpha',alphaValue)
                end
            else
                surf(X,Y,Z, 'FaceColor', color)
            end
        end
        if plotElementEdges
            plot3(X(:,1),Y(:,1),Z(:,1),'color',lineColor)
            plot3(X(:,end),Y(:,end),Z(:,end),'color',lineColor)

            plot3(X(1,:),Y(1,:),Z(1,:),'color',lineColor)
            plot3(X(end,:),Y(end,:),Z(end,:),'color',lineColor)
        end
        
    end
    for i = [1,Nxi]
        for j = [1,Neta]
            coeffs = reshape(c(1:3,i,j,[1,end]),3,2);
            plot3(coeffs(1,:),coeffs(2,:),coeffs(3,:),'color',lineColor)
        end
    end
% 
%     daspect([1 1 1]);
%     view(3)
%     camlight(110,30);
%     lighting phong;
%     set(gca, 'Color', 'none');
%     camlight('left')
    h = gcf;
end






