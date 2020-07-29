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

addpath ../../spheretri
nPoints = 100; 
tic;
[P,tri] = spheretri(nPoints);
P = 5*P;
noElems = size(tri,1)
% return
toc; 
trisurf(tri,P(:,1),P(:,2),P(:,3), 'FaceColor', getColor(1))
axis off
axis equal
camlight
TR = triangulation(tri,P);
% stlwrite(TR, ['../../../../../../../OneDrive/work/openfoam/BeTSSi/constant/triSurface/sphere_' num2str(noElems) '.stl']) 
 
