

scatteringCase = {'BI'}; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

model = {'IL'};

method = {'IE'};
% method = {'IENSG'};
% method = {'BEM'};
% BC = {'SHBC', 'SSBC','NNBC'};
% BC = {'SHBC'};
% BC = {'SSBC'};
BC = {'SSBC'};
% coreMethod = {'IGA'};
coreMethod = {'IGA'};
% coreMethod = {'hp_FEM'};
% {'IGA','hp_FEM','linear_FEM'}
c_f = 1524;
k_arr = 1;
omega_arr = k_arr*c_f;
f_arr = omega_arr/(2*pi); 


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

degreeElevArr = 0;
plotResultsInParaview = 0;
calculateFarFieldPattern = 0;
calculateVolumeError = 1;
calculateSurfaceError = 0;
plot2Dgeometry = 0;
plot3Dgeometry = 0;
scaleForErrorPlot = 0;
plotTS = 0;
plotMesh = 0;

collectIntoTasks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% coreMethod = {'IGA'};
% degreeElevArr = 0:2;
% 
% collectIntoTasks
