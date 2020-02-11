
if plot3Dgeometry
    figure('Color','white','name',['3D plot of geometry with mesh ' num2str(M)])
    if strcmp(method, 'IENSG')
        c_z = varCol{1}.c_z;
        c_xy = varCol{1}.c_xy;
        alignWithAxis = varCol{1}.alignWithAxis;
        x_0 = varCol{1}.x_0;
        ellipsoid = getEllipsoidalData(c_z,c_xy,c_xy,alignWithAxis, x_0);
        plotNURBS(ellipsoid,{'resolution',[20 40],'alphaValue',0.6,'color','blue'});
    end
    if strcmp(coreMethod, 'linear_FEM')
        resolution = [0,0,0];
    else 
        resolution = [40,40,40];
    end
    for j = 1:numel(varCol)
        fluid = varCol{j}.nurbs;
        for i = 1:numel(fluid)
            if strcmp(fluid{i}.type,'3Dvolume')
                if strcmp(coreMethod,'SEM')
                    varColPlot.plotAt = [0, 0;
                                         0, 0;
                                         1, 1];
                    plotLagrange(fluid{i},resolution, 1, 1.5*[44 77 32]/255, 0.8, NaN, varColPlot);
                else
                    plotNURBS(fluid{i},{'displayName',['Domain ' num2str(j) ', patch ' num2str(i)], ...
                                        'resolution',resolution,'alphaValue',0.8});
                end
            else
%                 plotNURBS(fluid{i},{'displayName',['Domain ' num2str(j) ', patch ' num2str(i)], 'resolution',resolution, ...
%                                                    'elementBasedSamples',true,'samplingDistance',0.1});
                plotNURBS(fluid{i},{'displayName',['Domain ' num2str(j) ', patch ' num2str(i)], 'resolution',resolution});
            end
        end
    end
    axis equal
    axis off
    set(gca, 'Color', 'none');
%     title(['Fluid 3D NURBS geometry. Mesh ' num2str(M)])
    drawnow
%     view(-70,30)
%     view(120,10)
    view(18,10) % BeTSSi
%     view(-72.3139,-1.5837)
%     view(18+180,10) % BeTSSi M1
%     view(11,25) % BeTSSi M2
%     view(62,43)
%     view(46,32) % BeTSSi M4
%     view(106,26) % sphere and cube
%     view(0,90) % above
    camproj('perspective')
%     camproj('orthographic')
    
    ax = gca;               % get the current axis
    ax.Clipping = 'off';    % turn clipping off

%     camlight

%     figureFullScreen(gcf)
% 	export_fig(['../../graphics/S1/SEM_p' num2str(degreeElev+4)], '-png', '-transparent', '-r200')
% 	export_fig(['../../graphics/S1/S13D_' coreMethod], '-png', '-transparent', '-r200')
% 	export_fig('../../graphics/TorusControlPts', '-png', '-transparent', '-r200')
% 	export_fig(['../../graphics/BCA/mesh' num2str(M)], '-png', '-transparent', '-r200')
% 	export_fig(['../../graphics/sphericalShell/S1patched'], '-png', '-transparent', '-r200')
% 	export_fig('C:\Users\jonvegar\Desktop\trond\mesh3D', '-png', '-transparent', '-r200')
% 	export_fig(['../../graphics/BCA/mesh' num2str(M)], '-png', '-transparent', '-r200')
% 	export_fig('../../graphics/M3/M3_above', '-png', '-transparent', '-r200')
% 	export_fig(['../../graphics/Cube_mesh' num2str(M)], '-png', '-transparent', '-r200')
% 	export_fig('../../graphics/M4', '-png', '-transparent', '-r200')
% 	export_fig('../../graphics/M2', '-png', '-transparent', '-r200')
% 	export_fig('../../graphics/M4_above', '-png', '-transparent', '-r200')
% 	export_fig('../../graphics/wineGlassFine', '-png', '-transparent', '-r200')
% 	export_fig('../../graphics/Shirron_above', '-png', '-transparent', '-r200')
% 	export_fig('../../graphics/M1_above', '-png', '-transparent', '-r200')
           

%     export_fig('../../graphics/sphericalShell/Sphere1controlPolygon' , '-png', '-transparent')
end   
if plot2Dgeometry && ~boundaryMethod
    fluid = varCol{1}.nurbs;
    solid = varCol{2}.nurbs;
    if useInnerFluidDomain
        fluid_i = varCol{3}.nurbs;
    end
%     figure(20+M)
    figure('Color','white','name',['3D plot of geometry with mesh ' num2str(M) ' - cross section'])
    hold on
    npts = 20;
    for i = 1:numel(fluid)
        switch model
%             case {'IL','S15'}
%                 if parm == 2 && i > 4
%                     continue
%                 end
%                 if parm == 1
%                     plotCuttingPlaneConstXi(fluid{i}, [npts npts],'xz',0);
%                 else
%                     plotCuttingPlaneConstEta(fluid{i}, [npts npts],'xz',0.5);
%                 end
            case {'IL','S15','SS','S1','S2','S3','PS','IL_P','SS_P','S1_P','S1_P2','S2_P','S3_P'}
                if parm == 2 && i > 4
                    continue
                end
                if parm == 1
                    plotCuttingPlaneConstXi(fluid{i}, [npts npts],'xz',0);
                else
                    plotCuttingPlaneConstEta(fluid{i}, [npts npts],'xz',0.5);
                end
            case {'BC','BC_P'}
    %             plotCuttingPlaneConstXi(fluid, [npts npts],'xz',0);
    %             plotCuttingPlaneConstXi(fluid, [npts npts],'xz',0.5);
                plotCuttingPlaneConstEta(fluid{i}, [npts npts],'yz',0.5);
            otherwise
                plotCuttingPlaneConstXi(fluid{i}, [npts npts],'xz',0.25);
    %             plotCuttingPlaneConstXi(fluid, [npts npts],'xz',0.75);
        end
        if useSolidDomain
            switch model
                case {'IL','S15','SS','S1','S2','S3','PS','IL_P','SS_P','S1_P','S1_P2','S2_P','S3_P'}
                    if parm == 1
                        plotCuttingPlaneConstXi(solid{i}, [npts npts],'xz',0, [0.8 0 0]);
                    else
                        plotCuttingPlaneConstEta(solid{i}, [npts npts],'xz',0.5, [0.8 0 0]);
                    end
                otherwise
                    plotCuttingPlaneConstXi(solid{i}, [npts npts],'xz',0.5, [0.8 0 0]);
            end
        end
        if useInnerFluidDomain
            switch model
                case {'IL','S15','SS','S1','S2','S3','PS','IL_P','SS_P','S1_P','S1_P2','S2_P','S3_P'}
                    if parm == 1
                        plotCuttingPlaneConstXi(fluid_i{i}, [npts npts],'xz',0);
                    else
                        plotCuttingPlaneConstEta(fluid_i{i}, [npts npts],'xz',0.5);
                    end
                otherwise
                    plotCuttingPlaneConstXi(fluid_i{i}, [npts npts],'xz',0.5);
            end
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
%     export_fig('../../graphics/sphericalShell/NNBC_FEMlinear6', '-pdf', '-transparent')
% 	export_fig(['../../graphics/MS_P/mesh' num2str(M) '_sphere'], '-pdf', '-transparent')
% 	export_fig(['../../graphics/MS_P/mesh' num2str(M)], '-pdf', '-transparent')
% 	export_fig(['../../graphics/BC/mesh' num2str(M) '_xi'], '-pdf', '-transparent')
% 	export_fig(['../../graphics/BC/mesh' num2str(M) '_eta'], '-pdf', '-transparent')
%     export_fig('../../graphics/sphericalShell/NNBC_IGAmesh4', '-pdf', '-transparent')
%     export_fig('../../graphics/M3/M3_above', '-png', '-transparent')
%     export_fig('../../graphics/M3/M3_M3_1', '-png', '-transparent')
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
        


