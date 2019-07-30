

scatteringCase = 'MS'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

model = 'M1'; % BeTSSi model 5A

f = 1e3; %[1e2 5e2 1e3];             % Frequency
omega = 2*pi*f;
k = omega/1500;

degree = 2;
alpha = (0:0.1:360)*pi/180;

plot3Dgeometry = 1;
plot2Dgeometry = 0;  % Plot cross section of mesh and geometry

solveForPtot = true;
plot3Dgeometry = 1;
method = {'BEM'};
formulation = {'CCBIE','CBM'};
formulation = {'CBM'};
M = 1;
loopParameters = {'formulation','M','method','f'};
collectIntoTasks
