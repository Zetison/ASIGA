

scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

model = 'IL';

method = 'IE';
formulation = 'BGU';
% method = {'IENSG'};
% method = {'BEM'};
% BC = {'SHBC', 'SSBC','NNBC'};
% BC = {'SHBC'};
% BC = {'SSBC'};
BC = 'NNBC';
% coreMethod = {'IGA'};
coreMethod = 'IGA';

c_f = 1524;
k = 2;
omega = k*c_f;
f = omega/(2*pi); 


M = 5; %3:6
% M = 5;
% N = 4;

alpha_s = pi;
beta_s = 0;
% alpha_s = 0;
% beta_s = -90*pi/180;
% alpha_s = 0*pi/180;
% beta_s = -90*pi/180;
% alpha_s = 180*pi/180;
% beta_s = 0*pi/180;

degreeElev = 0;
plotResultsInParaview = 0;
calculateFarFieldPattern = 0;
calculateVolumeError = 1;
calculateSurfaceError = 0;
plot2Dgeometry = 0;
plot3Dgeometry = 0;
scaleForErrorPlot = 0;
plotTS = 1;
plotMesh = 1;

collectIntoTasks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% coreMethod = {'IGA'};
% degreeElevArr = 0:2;
% 
% collectIntoTasks
