scatteringCase = 'MS';

model = 'Barrel';
method = {'BEM'};
formulation = 'CCBIE';

f = 1e3; % Frequency

M = 1;
parm = 1;
degree = 2;
alpha = (0:0.05:90)*pi/180;

loopParameters = {'M','parm','f','method'};
plot3Dgeometry = 1;
solveForPtot = true;
collectIntoTasks


method = {'KDT'};
solveForPtot = false;
formulation = 'MS1';
collectIntoTasks
