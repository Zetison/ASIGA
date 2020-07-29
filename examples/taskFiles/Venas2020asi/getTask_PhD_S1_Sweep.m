scatteringCase = 'Sweep';

model = 'S1';  % Spherical shell

% coreMethod = {'IGA', 'XI'};
coreMethod = {'IGA'};
varCol = setS1Parameters('double',1);
varCol{1}.meshFile = 'createNURBSmesh_EL';
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

k = 10.^linspace(-1,4,1000); % 10.^linspace(-1,2,1000)
% k = sort([eigenValues, k]);
% eigenValues = eigenValues*1500/(2*pi);
c_f = 1500;
f = k*c_f/(2*pi);

alpha = 0;
beta = pi;   

M = 1:3;
calculateSurfaceError = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RT simulation
coreMethod = {'IGA'};
method = {'RT'};
formulation = '';
M = 3;
prePlot.plot3Dgeometry = 0;
degree = 2;
calculateSurfaceError = 0;
computeCondNumber = false;
plotFarField = 1;
applyLoad = 'planeWave';
N = 3:6;
r = 1;
parm = 1;

loopParameters = {'method','M','parm','N','coreMethod'};
collectIntoTasks
