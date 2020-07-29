%% 3D Model Demo
% This is short demo that loads and renders a 3D model of a human femur. It
% showcases some of MATLAB's advanced graphics features, including lighting and
% specular reflectance.
close all
% Copyright 2011 The MathWorks, Inc.

% addpath ../export_fig

%% Load STL mesh
% Stereolithography (STL) files are a common format for storing mesh data. STL
% meshes are simply a collection of triangular faces. This type of model is very
% suitable for use with MATLAB's PATCH graphics object.

% Import an STL mesh, returning a PATCH-compatible face-vertex structure
if true
%     fv = stlread('../../../../../../../OneDrive/SINTEF/ATB/VTK_ATB_700_NEW/topo.stl');

    EdgeColor = [0,0,0];
%     EdgeColor = 'none';
    %% Render
    % The model is rendered with a PATCH graphics object. We also add some dynamic
    % lighting, and adjust the material properties to change the specular
    % highlighting.

    % patch(fv,'FaceColor',       [0.8 0.8 1.0], ...
    %          'EdgeColor',       'none',        ...
    %          'FaceLighting',    'gouraud',     ...
    %          'AmbientStrength', 0.15);
%     patch(fv,'FaceColor', [158,167,182]/255,'EdgeColor',EdgeColor);
    % patch(fv);

    % Add a camera light, and tone down the specular highlighting
    %     camproj('orthographic')
    % material('dull');

    % Fix the axes scaling, and set a nice view angle
    % figureFullScreen(gcf)
    % 
    % TR = triangulation(double(fv.faces),fv.vertices);
    % hold on  
    % P = incenter(TR);
    % F = faceNormal(TR);  
    % quiver3(P(:,1),P(:,2),P(:,3), ...
    %      F(:,1),F(:,2),F(:,3),0.5,'color','r');

    fv = stlread('../../../../../../../OneDrive/SINTEF/ATB/VTK_ATB_700_NEW/topo.stl');
    for i = 1:4
        figure(i)
        if i == 1
            patch(fv,'FaceColor', [102,120,70]/255,'EdgeColor',EdgeColor);
        end

        view(35.615966968232762,25.224770087692363) % BeTSSi
        ax = gca;               % get the current axis
        ax.Clipping = 'off';    % turn clipping off
        axis equal;
        axis off
        camproj('perspective')
        % camlight
        hold on
%         view([10.4767    8.4550         0   -9.4659;
%                -3.1984    5.0300    0.9046   -1.3681;
%                -6.7894   10.6774   -0.4262   95.4648;
%                      0         0         0    1.0000])
%         campos(1e6*[0.5960    6.6466    0.0003])
%         camtarget(1e6*[0.5958    6.6469    0.0001])
%         camup([0, 0, 1])
%         camva(11.3758)
        
        [fv.vertices,fv.faces] = refineTRI(fv.vertices,fv.faces);
        TR = triangulation(fv.faces,fv.vertices);
        name = ['../../../../../../../OneDrive/SINTEF/ATB/VTK_ATB_700_NEW/topo_' num2str(i)];
        format long
        stlwrite(TR, [name '.stl'],'text') 

%         options = struct('name',name,'celltype','VTK_TRIANGLE');
%         data.nodes = fv.vertices;
%         data.visElements = fv.faces;
%         makeVTKfile(data, options)
    end

    % TR = triangulation(double(fv.faces),fv.vertices);
    % hold on  
    % P = incenter(TR);
    % F = faceNormal(TR);  
    % quiver3(P(:,1),P(:,2),P(:,3), ...
    %      F(:,1),F(:,2),F(:,3),0.5,'color','r');
    %  
    %  
end
if false
    fv = stlread('../../../../../../../hugeFiles/VTK_ATB_700_NEW/bld.stl');

    % EdgeColor = [0,0,0];
    EdgeColor = 'none';
    %% Render
    % The model is rendered with a PATCH graphics object. We also add some dynamic
    % lighting, and adjust the material properties to change the specular
    % highlighting.

    % patch(fv,'FaceColor',       [0.8 0.8 1.0], ...
    %          'EdgeColor',       'none',        ...
    %          'FaceLighting',    'gouraud',     ...
    %          'AmbientStrength', 0.15);
    patch(fv,'FaceColor', [158,167,182]/255,'EdgeColor',EdgeColor);
    % patch(fv);

    % Add a camera light, and tone down the specular highlighting
    %     camproj('orthographic')
    % material('dull');

    % Fix the axes scaling, and set a nice view angle
    % figureFullScreen(gcf)
    % 
    % TR = triangulation(double(fv.faces),fv.vertices);
    % hold on  
    % P = incenter(TR);
    % F = faceNormal(TR);  
    % quiver3(P(:,1),P(:,2),P(:,3), ...
    %      F(:,1),F(:,2),F(:,3),0.5,'color','r');

    fv = stlread('../../../../../../../hugeFiles/VTK_ATB_700_NEW/topo.stl');

    patch(fv,'FaceColor', [102,120,70]/255,'EdgeColor',EdgeColor);

    view(35.615966968232762,25.224770087692363) % BeTSSi
    ax = gca;               % get the current axis
    ax.Clipping = 'off';    % turn clipping off
    axis equal;
    axis off
    camproj('perspective')
    % camlight


    % TR = triangulation(double(fv.faces),fv.vertices);
    % hold on  
    % P = incenter(TR);
    % F = faceNormal(TR);  
    % quiver3(P(:,1),P(:,2),P(:,3), ...
    %      F(:,1),F(:,2),F(:,3),0.5,'color','r');
    %  
    %  
end