
model = 'MS_P'; % Mock shell model after Ihlenburg

% method = {'IE','IENSG'};
% formulation = {'BGU','PGC','BGC'};
method = 'IE';
% formulation = {'BGU','PGU','PGC','BGC'};
% formulation = {'PGC','BGC'};
% formulation = {'BGC'};
formulation = {'BGU'};
% formulation = {'PGU','BGU'};
% formulation = {'PGU'};

IEbasis = {'Chebyshev', 'Lagrange', 'Bernstein'};
% IEbasis = {'Lagrange'};
% M = 3;
M = 1:6;
N = [1,2,3,6,9];
% N = 3;
degreeElevArr = 0;
plotResultsInParaview = 0;
calculateSurfaceError = 1;	% Only for spherical shell and if scatteringCase == 'Bi'
LpOrder = 2; % For error calculation in calcSurfError()

calculateVolumeError  = 1;	% Only for spherical shell and if scatteringCase == 'Bi'
calculateFarFieldPattern = 0;
computeCondNumber = 1;
applyLoad = 'radialPulsation';

k = 10;
c_f = 1500;

f = k*c_f/(2*pi);             % Frequency

alpha_s = 240*pi/180;
beta_s = 0*pi/180;     

plot2Dgeometry = 0;  % Plot cross section of mesh and geometry
plot3Dgeometry = 0;  % Plot cross section of mesh and geometry
loopParameters = {'M', 'N', 'formulation','IEbasis'};
% loopParameters = {'M'};
collectIntoTasks

method = {'BEM'};
formulation = {'CCBIE'};
% collectIntoTasks

%% Plot result in paraview
N = 6;
M = [1,3,4,6];
formulation = {'PGU'};
plot2Dgeometry = 1;  % Plot cross section of mesh and geometry
method = {'IE'};
plotResultsInParaview = 1;
calculateSurfaceError = 1;	% Only for spherical shell and if scatteringCase == 'Bi'
% collectIntoTasks
