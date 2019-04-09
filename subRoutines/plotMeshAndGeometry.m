
if plot3Dgeometry
    close all
    figure('Color','white','name',['3D plot of geometry with mesh ' num2str(M)])
    for i = 1:numel(fluid)
        if boundaryMethod
            switch method
                case 'IENSG'
                    c_z = varCol.c_z;
                    c_xy = varCol.c_xy;
                    alignWithAxis = varCol.alignWithAxis;
                    x_0 = varCol.x_0;
                    ellipsoid = getEllipsoidalData(c_z,c_xy,c_xy,alignWithAxis, x_0);
                    plotNURBS(ellipsoid,{'resolution',[20 40],'alphaValue',0.6,'color','blue'});

                    hold on
                    plotNURBS(fluid{i},{'resolution',[40 40],'alphaValue',0.8});
                    hold off
                case {'BEM', 'KDT','MFS','RT','BA'}
                    if strcmp(coreMethod, 'linear_FEM')
                        resolution = [0,0];
                    else 
                        resolution = [10,10];
                    end
                    if 0
                        colorFun = @(v) abs(varCol.analytic(v));
%                         colorFun = @(v) real(varCol.analytic(v));
                        plotNURBS(fluid{i},{'resolution',[10 10 10],'colorFun',colorFun,'elementBasedSamples',true, ...
                                                        'samplingDistance',0.1});
                    else
                        plotNURBS(fluid{i},{'resolution',[100 100 10]});
                    end
            end
        else
            varColPlot.plotAt = [0, 0;
                                 0, 0;
                                 1, 1];
            if strcmp(coreMethod,'SEM')
%                     colorFun = @(x) log10(abs(norm2(x)));
                plotLagrange(fluid{i},[1000 1000 1000], 1, 1.5*[44 77 32]/255, 0.8, NaN, varColPlot);
            else
                plotAt = [0, 0;
                         0, 0;
                         0, 1];
                if boundaryMethod
                    plotNURBS(fluid{i},{'resolution',[10 10 10]});
                else
                    colorFun = @(x) log10(abs(norm(x)-2)/2);
                    if false
                        plotNURBS(fluid{i},{'resolution',[100 100 100],'alphaValue',0.8,'colorFun',colorFun,'plotAt',plotAt});
                        caxis([-17,1])
                        plotControlPts(fluid{i})
                    else
                        plotNURBS(fluid{i},{'resolution',[10 10 10],'alphaValue',1,'plotAt',plotAt});
                    end
                end
            end
%                 plotNURBS(fluid{i},[100 100 0], 1, 1.5*[44 77 32]/255, 1);

            if 0 % test artificial boundary
%             nurbs = extractOuterSurface(fluid);
%             x_0 = varCol.x_0;
%             c_z = 34; % 30
%             c_xy = c_z/5; % 2.5, 3.75
%             a = c_z;
%             b = c_xy;
%             c = c_xy;
%             varCol.colorFun = @(v) abs((v(1)-x_0(1))^2/a^2 + (v(2)-x_0(2))^2/b^2 + (v(3)-x_0(3))^2/c^2 - 1);
%             npts = 100;
%             plotNURBS(nurbs,[npts npts], 1, 1.5*[44 77 32]/255, 1, NaN, varCol);
            end
        end
    end
    axis equal
    axis off
    set(gca, 'Color', 'none');
%     title(['Fluid 3D NURBS geometry. Mesh ' num2str(M)])
    drawnow
    varCol.dimension = 1;
%     view(-70,30)
%     view(120,10)
    view(18,10)
%     view(0,0)
    camproj('perspective')
    ax = gca;               % get the current axis
    ax.Clipping = 'off';    % turn clipping off
    camlight
%     keyboard 
%     plotControlPts(fluid,'red','red','red','red')
    figureFullScreen(gcf)
% 	export_fig(['../graphics/S1/SEM_p' num2str(degreeElev+4)], '-png', '-transparent', '-r300')
% 	export_fig(['../graphics/S1/S13D_' coreMethod], '-png', '-transparent', '-r300')
% 	export_fig(['../graphics/BCA/mesh' num2str(M)], '-png', '-transparent', '-r300')
% 	export_fig('C:\Users\jonvegar\Desktop\trond\mesh3D', '-png', '-transparent', '-r300')
           

%     export_fig(['../graphics/sphericalShell/patchedShellMesh' num2str(M)], '-png', '-transparent')
end   
if plot2Dgeometry && ~boundaryMethod
%     figure(20+M)
    figure('Color','white','name',['3D plot of geometry with mesh ' num2str(M) ' - cross section'])
    hold on
    npts = 20;
    for i = 1:numel(fluid)
        switch model
            case {'IL','SS','S1','S2','S3','PS','IL_P','SS_P','S1_P','S1_P2','S2_P','S3_P'}
                if parm == 2 && i > 4
                    continue
                end
                plotCuttingPlaneConstEta(fluid{i}, [npts npts],'xz',0.5);
            case {'BC','BC_P'}
    %             plotCuttingPlaneConstXi(fluid, [npts npts],'xz',0);
    %             plotCuttingPlaneConstXi(fluid, [npts npts],'xz',0.5);
                plotCuttingPlaneConstEta(fluid{i}, [npts npts],'yz',0.5);
            otherwise
                plotCuttingPlaneConstXi(fluid{i}, [npts npts],'xz',0.25);
    %             plotCuttingPlaneConstXi(fluid, [npts npts],'xz',0.75);
        end
        if useInnerFluidDomain
            plotCuttingPlaneConstXi(fluid_i{1}, [npts npts],'xz',0.25);
%             plotCuttingPlaneConstXi(fluid_i{1}, [npts npts],'xz',0.5);
    %         plotCuttingPlaneConstEta(fluid_i, [npts npts],'yz',0.5);
        end
        if useSolidDomain
            plotCuttingPlaneConstXi(solid{1}, [npts npts],'xz',0.25, [0.8 0 0]);
%             plotCuttingPlaneConstXi(solid{1}, [npts npts],'xz',0.5, [0.8 0 0]);
    %         plotCuttingPlaneConstEta(solid, [npts npts],'yz',0.5, [0.8 0 0]);
        end
    end
    hold off
    axis equal
    axis off
    set(gca, 'Color', 'none');
    set(gca, 'Color', 'white');
    set(gcf,'color','w');
%     title(['Cross section of Fluid 3D NURBS geometry. Mesh ' num2str(i_mesh)])
    drawnow
%     export_fig('../graphics/sphericalShell/NNBC_FEMlinear6', '-pdf', '-transparent')
% 	export_fig(['../graphics/MS_P/mesh' num2str(M) '_sphere'], '-pdf', '-transparent')
% 	export_fig(['../graphics/MS_P/mesh' num2str(M)], '-pdf', '-transparent')
% 	export_fig(['../graphics/BC/mesh' num2str(M) '_xi'], '-pdf', '-transparent')
% 	export_fig(['../graphics/BC/mesh' num2str(M) '_eta'], '-pdf', '-transparent')
%     export_fig('../graphics/sphericalShell/NNBC_IGAmesh4', '-pdf', '-transparent')
%     export_fig('../graphics/M3/M3_above', '-png', '-transparent')
%     export_fig('../graphics/M3/M3_M3_1', '-png', '-transparent')
end
% keyboard




return
Eps = 4e-2;
for i = 1:size(varCol.controlPts,1)
    P = varCol.controlPts(i,:);
    r = norm(P);
    x = P(1);
    y = P(2);
    z = P(3);
    theta = acos(z/r);
    phi = atan2(y,x);
    T = r*[cos(theta)*cos(phi),cos(theta)*sin(phi),-sin(theta)];
    T2 = r*[-sin(theta)*sin(phi),sin(theta)*cos(phi),0];
%     Dr = Eps*T*rand() + Eps*T2*rand();
    Dr = (r+Eps*i)*[sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)];
    text(Dr(1),Dr(2),Dr(3),num2str(i),'color','blue')
%     text(P(1),P(2),P(3),num2str(i),'color','red')
end
p = 2;
for i = 1:2:size(fluid.coeffs,2)-2
    for j = 1:2:size(fluid.coeffs,3)-2
        pts = reshape(fluid.coeffs(1:3,i:i+p,j:j+p),3,(p+1)^2);
        maxX = max(pts,[],2);
        minX = min(pts,[],2);
        X = [minX(1), maxX(1);
             minX(1), maxX(1)];
        Y = [minX(2), minX(2);
             maxX(2), maxX(2)];
        Z = [minX(3), minX(3);
             minX(3), minX(3)];
        surf(X,Y,Z, 'FaceColor', 'red','FaceAlpha',1, 'FaceLighting', 'none','FaceAlpha',0.8)
        
        X = [minX(1), maxX(1);
             minX(1), maxX(1)];
        Y = [minX(2), minX(2);
             maxX(2), maxX(2)];
        Z = [maxX(3), maxX(3);
             maxX(3), maxX(3)];
        surf(X,Y,Z, 'FaceColor', 'red','FaceAlpha',1, 'FaceLighting', 'none','FaceAlpha',0.8)
        
        X = [minX(1), maxX(1);
             minX(1), maxX(1)];
        Z = [minX(3), minX(3);
             maxX(3), maxX(3)];
        Y = [minX(2), minX(2);
             minX(2), minX(2)];
        surf(X,Y,Z, 'FaceColor', 'red','FaceAlpha',1, 'FaceLighting', 'none','FaceAlpha',0.8)
        
        X = [minX(1), maxX(1);
             minX(1), maxX(1)];
        Z = [minX(3), minX(3);
             maxX(3), maxX(3)];
        Y = [maxX(2), maxX(2);
             maxX(2), maxX(2)];
        surf(X,Y,Z, 'FaceColor', 'red','FaceAlpha',1, 'FaceLighting', 'none','FaceAlpha',0.8)
        
        Z = [minX(3), maxX(3);
             minX(3), maxX(3)];
        Y = [minX(2), minX(2);
             maxX(2), maxX(2)];
        X = [minX(1), minX(1);
             minX(1), minX(1)];
        surf(X,Y,Z, 'FaceColor', 'red','FaceAlpha',1, 'FaceLighting', 'none','FaceAlpha',0.8)
        
        Z = [minX(3), maxX(3);
             minX(3), maxX(3)];
        Y = [minX(2), minX(2);
             maxX(2), maxX(2)];
        X = [maxX(1), maxX(1);
             maxX(1), maxX(1)];
        surf(X,Y,Z, 'FaceColor', 'red','FaceAlpha',1, 'FaceLighting', 'none','FaceAlpha',0.8)
        drawnow
    end
end
        


