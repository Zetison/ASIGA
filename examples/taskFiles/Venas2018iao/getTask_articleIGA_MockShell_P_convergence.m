%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This study correspond to Figure 21-23 in Venas2018iao
% Venas2018iao is available at https://doi.org/10.1016/j.cma.2018.02.015 (open access version at http://hdl.handle.net/11250/2493754)


model = 'MS_P'; % Mock shell model after Ihlenburg

% method = {'IE','IENSG'};
% formulation = {'BGU','PGC','BGC'};
method = 'IE';
formulation = {'BGU','PGU','PGC','BGC'};
% formulation = {'BGU'};

IEbasis = {'Chebyshev','Bernstein','Lagrange'};
% IEbasis = {'Bernstein'};
M = 1:6;
M = 1:2;
N = [1,2,3,6,9];
% N = 1;
degree = 2;
runTasksInParallel = 1;
plotResultsInParaview = 0;
calculateSurfaceError = 1;	% Only for spherical shell and if scatteringCase == 'Bi'
LpOrder = 2; % For error calculation in calcSurfError()

calculateVolumeError  = 1;	% Only for spherical shell and if scatteringCase == 'Bi'
calculateFarFieldPattern = 0;
computeCondNumber = 1;
applyLoad = 'pointPulsation';
BC = 'NBC';

k = 10;
c_f = 1500;

f = k*c_f/(2*pi);             % Frequency

prePlot.plot2Dgeometry = 0;  % Plot cross section of mesh and geometry
prePlot.plot3Dgeometry = 0;  % Plot cross section of mesh and geometry
loopParameters = {'M', 'N', 'formulation','IEbasis'};
% loopParameters = {'M'};
collectIntoTasks
