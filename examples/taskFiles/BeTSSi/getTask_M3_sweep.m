

misc.scatteringCase = 'Sweep'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

misc.model = 'M3'; % BeTSSi misc.model 5A

misc.method = {'IE','IENSG'};
% misc.method = {'IENSG'};
formulation = 'BGU';

f = linspace(1e2,2e2,2); %[1e2 5e2 1e3];             % Frequency
omega = 2*pi*f;
k = omega/1500;
N = 3;

M = 1;
degree = 2:3;
loopParameters = {'misc.f','M','degree','N','alpha','misc.method'};
alpha = [90,180]*pi/180;

alpha_s = 240*pi/180;
beta_s = 0*pi/180;        
% beta_s = 30*pi/180;
prePlot.plot3Dgeometry = 1;

prePlot.plot2Dgeometry = 1;  % Plot cross section of mesh and geometry
% collectIntoTasks

misc.method = {'KDT'};
collectIntoTasks
