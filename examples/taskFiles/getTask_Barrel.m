scatteringCase = 'MS';

model = 'Barrel';
method = {'BEM'};
formulation = {'GBM','GCBIE','CBM','CCBIE'};
BC = 'SHC';

f = 1.5e3; % Frequency

parms = setBarrelParameters();

M = 1;
parm = 2;
degree = 2;
alpha = (0:0.05:180)*pi/180;

loopParameters = {'M','parm','f','method','formulation'};
prePlot.plot3Dgeometry = false;
solveForPtot = true;
collectIntoTasks


method = {'KDT'};
solveForPtot = false;
formulation = {'MS1'};
% collectIntoTasks
