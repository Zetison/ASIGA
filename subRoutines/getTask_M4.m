

scatteringCase = 'MS'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

model = 'M4'; % BeTSSi model 5A and BeTSSi model 5B

method = {'BEM'};
formulation = {'CBM','CCBIE'};
% formulation = {'CCBIE'};

f = 10e3;             % Frequency
f = 1e3;             % Frequency

M = 3:6;
% M = 1;
parm = [2,1];
% parm = 2;
degree = 4;
beta = 30*pi/180;
beta = 0;
alpha = (0:0.05:360)*pi/180;

loopParameters = {'M','parm','f','method','formulation'};
plot3Dgeometry = 0;
solveForPtot = true;
% collectIntoTasks


method = {'KDT'};
M = 3;
parm = 2;
solveForPtot = false;
extraGP = 0; % extra quadrature points
extraGPBEM = 10; % extra quadrature points around singularities for BEM formulations
agpBEM = 0.6; % parameter for adaptiv Gauss point integration around singularities for BEM formulations

formulation = {'MS1'};
collectIntoTasks
