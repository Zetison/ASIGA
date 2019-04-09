

scatteringCase = {'BI'}; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

model = {'S1_P'};

method = {'IE','IENSG'};
% method = {'IENSG'};
formulation = {'PGU','BGU','PGC','BGC'};
BC = {'NBC'};
plot3Dgeometry = 0;
plot2Dgeometry = 0;
calculateFarFieldPattern = false;
applyLoad = 'radialPulsation';
k = 10;
c_f = 1500;

f_arr = k*c_f/(2*pi);             % Frequency

M = 1:6; %3:6
N = [1,2,3,6,9];

alpha_s = 180*pi/180;
beta_s = 0*pi/180;


degreeElevArr = 0;

calculateSurfaceError = 1;
calculateVolumeError  = 1;
LpOrder = Inf;
collectIntoTasks

%% Plot result in paraview
% N = 6;
% M = 5;
% formulation = {'PGU'};
% method = {'IE'};
% plotResultsInParaview = 1;
% calculateSurfaceError = 1;	% Only for spherical shell and if scatteringCase == 'Bi'
% collectIntoTasks
