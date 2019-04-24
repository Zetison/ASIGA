

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
% formulation = {'PGU', 'PGC'};
formulation = {'PGU'};
coreMethod = {'IGA'};
if strcmp(method, 'BEM')
    formulation = {'CCBIE'};
%     formulation = {'CHBIE'};
%     formulation = {'CBM'};
end

f_arr = 1.5e3;             % Frequency

M = 5:6; %3:6
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
calculateFarFieldPattern = 1;
calculateVolumeError = 1;
calculateSurfaceError = 0;
plot2Dgeometry = 1;
plot3Dgeometry = 0;

collectIntoTasks
