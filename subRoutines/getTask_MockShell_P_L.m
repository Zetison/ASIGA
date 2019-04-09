
model = 'MS_P'; % Mock shell model after Ihlenburg

% method = {'IE','IENSG'};
% formulation = {'BGU','PGC','BGC'};
method = 'IE';
% formulation = {'BGU','PGU','PGC','BGC'};
% formulation = {'PGC','BGC'};
% formulation = {'BGC'};
formulation = 'BGU';
% formulation = {'PGU','BGU'};
% formulation = {'PGU'};

IEbasis = 'Chebyshev';
% IEbasis = {'Lagrange'};
% M = 1:6;
M = 6;
N = 3;
% N = 1:2:5;
degreeElevArr = 0;
plotResultsInParaview = 1;
plotMesh = 1;
calculateSurfaceError = 1;	% Only for spherical shell and if scatteringCase == 'Bi'
LpOrder = Inf; % For error calculation in calcSurfError()

calculateVolumeError  = 1;	% Only for spherical shell and if scatteringCase == 'Bi'
calculateFarFieldPattern = 0;
computeCondNumber = 0;
applyLoad = 'radialPulsation';
BC = 'NBC';

% parm = [1,4,7,10];
% BC_length = 5.075*(6.5+6.5+42+7.0)/3.5 - 2*5.075;
% parm = BC_length; % *(6.5+6.5+42+7.0)/4
parm = 1;

k = 10;
c_f = 1500;

f = k*c_f/(2*pi);             % Frequency
f = 5e3;

alpha_s = 240*pi/180;
beta_s = 0*pi/180;     

plot2Dgeometry = 0;  % Plot cross section of mesh and geometry
plot3Dgeometry = 0;  % Plot cross section of mesh and geometry
loopParameters = {'parm', 'M', 'N'};
% loopParameters = {'M'};
collectIntoTasks
