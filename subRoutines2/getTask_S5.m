scatteringCase = 'Sweep';

model = 'S5';  % Spherical shell

% coreMethod = {'IGA', 'XI'};
coreMethod = {'IGA'};
method = {'BEM'};
if strcmp(method, 'BEM')
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

k = 10.^linspace(-1,1,3000); % 10.^linspace(-1,2,1000)
k = sort([eigenValues, k]);
eigenValues = eigenValues*1500/(2*pi);
c_f = 1500;
f = k*c_f/(2*pi);


alpha_s = 0;
beta_s = pi;  

alpha = alpha_s;
beta = beta_s;   

M = 1:3;
calculateSurfaceError = 1;

% collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MFS simulation
method = {'MFS'};
coreMethod = 'IGA';
formulation = '';
M = 5;
degreeElev = 0;
calculateSurfaceError = 0;
computeCondNumber = false;
loopParameters = {'f','M','method'};
% collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% KDT simulation
method = {'KDT'};
formulation = '';
M = 3;
plot3Dgeometry = 0;
degreeElev = 0;
calculateSurfaceError = 0;
computeCondNumber = false;
loopParameters = {'M','method'};
collectIntoTasks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RT simulation
method = {'RT'};
formulation = '';
M = 3;
plot3Dgeometry = 0;
degreeElev = 0;
calculateSurfaceError = 0;
computeCondNumber = false;
plotFarField = 1;
applyLoad = 'planeWave';
parm = 3:6;
r = 2;

loopParameters = {'parm','M','method'};
% collectIntoTasks
