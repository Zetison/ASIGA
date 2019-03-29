

scatteringCase = {'BI'}; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

model = {'S2'};

method = {'IE'};
% method = {'IENSG'};
% method = {'BEM'};
% BC = {'SHBC', 'SSBC','NNBC'};
BC = {'SHBC'};
% BC = {'SSBC'};
% BC = {'NNBC'};
% coreMethod = {'IGA'};
coreMethod = {'IGA', 'C0_IGA', 'hp_FEM','h_FEM','linear_FEM'};
if strcmp(method, 'BEM')
    formulation = {'CCBIE'};
%     formulation = {'CHBIE'};
%     formulation = {'CBM'};
end

f_arr = 1e3;             % Frequency

M = 1:7; %3:6
% M = 5;
% N = 4;

alpha_s = 240*pi/180;
beta_s = 30*pi/180;
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

collectIntoTasks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
coreMethod = {'IGA'};
degreeElevArr = 0:2;

collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% method = {'BEM'};
% coreMethod = {'IGA', 'C0_IGA', 'hp_FEM'};
% formulation = {'CCBIE'};
% M = 2:5;
% useSolidDomain = false;
% degreeElevArr = 0;
% calculateVolumeError = 0;
% calculateSurfaceError = 1;
% 
% collectIntoTasks
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% coreMethod = {'IGA'};
% degreeElevArr = 1:4;
% 
% collectIntoTasks