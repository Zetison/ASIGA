misc.scatteringCase = 'Sweep';

misc.model = 'S1';  % Spherical shell

% misc.coreMethod = {'IGA', 'XI'};
misc.coreMethod = {'IGA'};
misc.method = {'BEM'};
if strcmp(misc.method, 'BEM')
    formulation = {'CCBIE', 'CBM', 'CHBIE'};
end


f = linspace(0.05,2e3,50);             % Wave number for outer fluid domain
% Eigenfrequencies from Zheng2015itb (note that they are multiplied by two, to
% compensate for the half size of the sphere)
eigenValues = [3.141592653589794 % pi
                   6.283185307179587 % 2*pi
                   9.424777960769379 % 3*pi
                   4.493409457909065
                   7.725251836937707
                   5.763459196894550
                   9.095011330476353
                   6.987932000500519
                   8.182561452571243
                   9.355812111042747]';

k = 10.^linspace(-1,4,1000); % 10.^linspace(-1,2,1000)
% k = sort([eigenValues, k]);
% eigenValues = eigenValues*1500/(2*pi);
c_f = 1500;
f = k*c_f/(2*pi);

alpha = 0;
beta = pi;   

M = 1:3;
err.calculateSurfaceError = 1;

% collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MFS simulation
misc.method = {'MFS'};
misc.coreMethod = 'IGA';
formulation = '';
M = 5;
degree = 2;
err.calculateSurfaceError = 0;
computeCondNumber = false;
loopParameters = {'f','M','misc.method'};
% collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% KDT simulation
misc.method = {'KDT'};
formulation = '';
prePlot.plot3Dgeometry = 0;
degree = 4;
err.calculateSurfaceError = 0;
computeCondNumber = false;
loopParameters = {'misc.method','M','parm','misc.coreMethod'};
misc.coreMethod = {'linear_FEM'};
% misc.coreMethod = {'IGA'};
M = [6,8,10];
M = 10;
parm = 2;
% k = 10.^linspace(-1,4,1000);
k = 1000;
c_f = 1500;
f = k*c_f/(2*pi);
collectIntoTasks


% misc.coreMethod = {'IGA'};
% % M = 3;
% collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RT simulation
misc.method = {'RT'};
formulation = '';
M = 3;
prePlot.plot3Dgeometry = 0;
degree = 2;
err.calculateSurfaceError = 0;
computeCondNumber = false;
plotFarField = 1;
misc.applyLoad = 'planeWave';
N = 3:6;
r = 1;
parm = 1;

loopParameters = {'misc.method','M','parm','N','misc.coreMethod'};
% collectIntoTasks
