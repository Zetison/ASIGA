

scatteringCase = 'MS'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

model = 'M4'; % BeTSSi model 5A and BeTSSi model 5B

method = {'BEM'};
formulation = {'CBM','CCBIE'};
% formulation = {'CCBIE'};

f = 10e3;             % Frequency

M = 3:6;
M = 7;
parm = [2,1];
% parm = 2;
degree = 4;
beta = 30*pi/180;
alpha = (0:0.05:360)*pi/180;

loopParameters = {'M','parm','f','method','formulation'};
plot3Dgeometry = 0;
solveForPtot = true;
collectIntoTasks


method = {'KDT'};
solveForPtot = false;
formulation = {''};
collectIntoTasks
