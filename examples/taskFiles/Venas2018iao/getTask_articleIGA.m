

misc.scatteringCase = {'BI'}; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

misc.model = {'S2'};

misc.method = {'IE'};
% misc.method = {'IENSG'};
% misc.method = {'BEM'};
% BC = {'SHBC', 'SSBC','NNBC'};
BC = {'SHBC'};
% BC = {'SSBC'};
% BC = {'NNBC'};
% misc.coreMethod = {'IGA'};
misc.coreMethod = {'IGA', 'C0_IGA', 'hp_FEM','h_FEM','linear_FEM'};
if strcmp(misc.method, 'BEM')
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
err.calculateSurfaceError = 0;
plot2Dgeometry = 0;
plot3Dgeometry = 0;

collectIntoTasks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
misc.coreMethod = {'IGA'};
degreeElevArr = 0:2;

collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% misc.method = {'BEM'};
% misc.coreMethod = {'IGA', 'C0_IGA', 'hp_FEM'};
% formulation = {'CCBIE'};
% M = 2:5;
% useSolidDomain = false;
% degreeElevArr = 0;
% calculateVolumeError = 0;
% err.calculateSurfaceError = 1;
% 
% collectIntoTasks
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% misc.coreMethod = {'IGA'};
% degreeElevArr = 1:4;
% 
% collectIntoTasks