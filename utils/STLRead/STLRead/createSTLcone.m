%% 3D Model Demo
% This is short demo that loads and renders a 3D model of a human femur. It
% showcases some of MATLAB's advanced graphics features, including lighting and
% specular reflectance.
close all
clear all
% Copyright 2011 The MathWorks, Inc.

% addpath ../export_fig

%% Load STL mesh
% Stereolithography (STL) files are a common format for storing mesh data. STL
% meshes are simply a collection of triangular faces. This type of model is very
% suitable for use with MATLAB's PATCH graphics object.

% Import an STL mesh, returning a PATCH-compatible face-vertex structure
axisCenter = [595866, 6646920, 116.5];
% axisCenter = [0,0,0.5];
angleOfAttack = 40*pi/180;
h_max = 500;
% h_max = 1;
n1 = 300;
n2 = 100;
theta = linspace(0,2*pi,n1).';
% theta = theta(1:end-1);
h = linspace(axisCenter(3),h_max,n2)-axisCenter(3);
r = h/tan(angleOfAttack);
[H,THETA] = meshgrid(h,theta);
H = reshape(H,n1*n2,1);
THETA = reshape(THETA,n1*n2,1);
X = axisCenter(1)+cos(theta)*r;
Y = axisCenter(2)+sin(theta)*r;
Z = ones(size(theta))*(h+axisCenter(3));
X = reshape(X,n1*n2,1);
Y = reshape(Y,n1*n2,1);
Z = reshape(Z,n1*n2,1);
% plot3(X,Y,Z,'*')
hold on
X = [X,Y,Z];
DT = delaunay([H,THETA]);
TR = triangulation(DT,X);
EdgeColor = 'black';
trimesh(TR.ConnectivityList,X(:,1),X(:,2),X(:,3),'FaceColor', [102,120,70]/255,'EdgeColor',EdgeColor)
% triplot(TR.ConnectivityList,H,THETA)
% return
view(35.615966968232762,45.224770087692363) % BeTSSi
ax = gca;               % get the current axis
ax.Clipping = 'off';    % turn clipping off
axis equal;
axis off
camproj('perspective')
camlight
drawnow
stlwrite(TR, '../../../../../../../hugeFiles/VTK_ATB_700_NEW/cone.stl') 
 
 