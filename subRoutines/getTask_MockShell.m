

scatteringCase = 'BI'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

model = 'MS'; % Mock shell model after Ihlenburg

calculateVolumeError  = 0;	% Only for spherical shell and if scatteringCase == 'Bi'
calculateFarFieldPattern = 0;

k = 10;
c_f = 1500;

f = k*c_f/(2*pi);             % Frequency

alpha_s = 240*pi/180;
beta_s = 0*pi/180;     

plot2Dgeometry = 0;  % Plot cross section of mesh and geometry
plot3Dgeometry = 0;  % Plot cross section of mesh and geometry

%% Plot result in paraview
N = 4;
M = 5;
formulation = 'BGU';
method = 'IE';
BC = 'SSBC';

plotResultsInParaview = 1;
calculateSurfaceError = 0;	% Only for spherical shell and if scatteringCase == 'Bi'


loopParameters = {'M'};
collectIntoTasks
