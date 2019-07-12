

scatteringCase = 'MS'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

model = 'M4'; % BeTSSi model 5A and BeTSSi model 5B

method = {'BEM'};
formulation = {'CBM','CCBIE'};

f = [1e3, 3e3, 10e3];             % Frequency
% f = 1e3;             % Frequency

M = 3:6;
% M = 3:4;
% parm = 2;
parm = [2,1];
degree = 2;
beta_s = 30*pi/180;
alpha = (0:0.05:360)*pi/180;

loopParameters = {'M','parm','f','method','formulation'};
plot3Dgeometry = 0;
solveForPtot = true;
collectIntoTasks


method = {'KDT'};
solveForPtot = false;
formulation = {''};
collectIntoTasks
