% Copyright 2011 The MathWorks, Inc.
clear all
close all

%% Load STL mesh
% Stereolithography (STL) files are a common format for storing mesh data. STL
% meshes are simply a collection of triangular faces. This type of model is very
% suitable for use with MATLAB's PATCH graphics object.

% Import an STL mesh, returning a PATCH-compatible face-vertex structure
fv = stlread('../../../../rhinoceros/BeTSSi/BeTSSi_mod_10cm.stl');
% fv = stlread('../../../../rhinoceros/test.stl');
% fv = stlread('../../../../../FFI/BeTSSiIIb/Model_Files/outer_hull/hull_20cmMesh.stl');


%% Render
% The model is rendered with a PATCH graphics object. We also add some dynamic
% lighting, and adjust the material properties to change the specular
% highlighting.

patch(fv,'FaceColor',       [0.8 0.8 1.0], ...
         'FaceLighting',    'gouraud',     ...
         'AmbientStrength', 0.15);
%          'EdgeColor',       'none',        ...

% Add a camera light, and tone down the specular highlighting
camlight('headlight');
material('dull');
ax = gca;               % get the current axis
ax.Clipping = 'off';    % turn clipping off
axis equal
axis off

% Fix the axes scaling, and set a nice view angle
% axis('image');
view([-135 35]);

