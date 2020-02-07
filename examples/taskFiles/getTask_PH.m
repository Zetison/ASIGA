

scatteringCase = 'MS'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

model = 'PH'; % BeTSSi model 5A

% method = {'IE'};
method = {'IENSG','IE'};
formulation = 'BGU';
BC = 'SHBC';

f = 2e2; %[1e2 5e2 1e3];             % Frequency

M = 1;
loopParameters = {'M','method'};
alpha = (0:0.05:360)*pi/180;
degree = 2;

alpha_s = 240*pi/180;
beta_s = 0*pi/180;        
% beta_s = 30*pi/180;
plot3Dgeometry = 0;

plot2Dgeometry = 0;  % Plot cross section of mesh and geometry
% collectIntoTasks

method = {'KDT'};
M = 3:4;
% M = 1;
collectIntoTasks